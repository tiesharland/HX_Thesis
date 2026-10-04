function results = run_lumped_model(e, r, hx_theta, fan, fpr_init, ...
    M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
    d2_init, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, ...
    solve_T, save_results)
% run_lumped_model  Forward-pass pressure/mass-flow solve for the lumped
% (non-discretised) HX model.
%
% solve_T: bookkeeping only for now (not yet implemented), kept for
% parity with run_disc_model_fwdpass -- see that function's header.
%
% save_results: if true, saves this run's results to disk under
% ./results, keyed into the SAME "physical input configuration" index
% table (results_index.mat) used by run_disc_model_fwdpass, so the same
% id can refer to the lumped model's results and the discretised model's
% results (1D/2D, DNS or not) at the same physical inputs. The model
% variant is encoded in the results filename only: where
% run_disc_model_fwdpass uses "a<N_segments>[-c<N_cool_seg>][_dns]_<id>.mat",
% this function always uses "lumped_<id>.mat" in its place, since there is
% no N_segments/N_cool_seg/use_DNS choice here.
%
% IMPORTANT: exactly as in run_disc_model_fwdpass, this function always
% runs the full solve -- it never skips computation, even if a matching
% results file already exists. The matching-against-results_index step is
% purely to decide which id (and therefore which results file) this run
% should be saved as/into.
%
% overwrite (hard-coded local flag, just below): same meaning as in
% run_disc_model_fwdpass.
%   true  (default) -> overwrite lumped_<id>.mat in place.
%   false            -> keep the existing file(s) and version instead,
%                        e.g. lumped_3.mat, then
%                        lumped_3-1.mat, lumped_3-2.mat, ...
%
% The table's columns are every physical input this function takes
% (fpr_init is recorded as the caller's input value, not the possibly
% fan-adjusted value computed below), plus solve_T.

overwrite = true;   % hard-coded; set to false to version instead of overwrite

results_dir = 'results';
index_file  = fullfile(results_dir, 'results_index.mat');

fpr_init_arg = fpr_init;   % capture the input BEFORE it may be reassigned
                           % below, so the saved "physical input" key
                           % reflects what was passed in, not the
                           % fan-adjusted value.

t_compute_start = tic;

counter = 1;
x=1;

%%%--- SHARED: FREESTREAM ---%%%

[P_inf, T_inf, Rho_inf, Gamma_inf, Cp_inf, mu_inf, P_inf_tot, M_inf, a_inf] = ...
    Freestream(h, p11, R, t11, V_inf);

P_inf_tot = P_inf + (0.5*Rho_inf*V_inf*V_inf);
a_inf     = sqrt(Gamma_inf*R*T_inf);
M_inf     = V_inf/a_inf;

%%%--- SHARED: PROPELLER ---%%%

[J_val, Cp_val, eff_val, D, J, Cp, V_ind_prop, dp_tot_prop, T_net, dia_prop, ...
    M1, P1, P1_0, T1, T1_0, Rho_1, v1, Re_1] = ...
    propeller(flight_phase, V_inf, Rho_inf, P_inf, a_inf, T_inf, mu_inf);

%%%--- SHARED: DIFFUSER ---%%%

A2_init = pi*(d2_init^2)/4; %m^2
% AR_init = AR_diff_arr; %initial assumption for area ratio (A2/A1)


A2             = pi*d2_init*d2_init/4;
A_frontal_target = A2*AR_diff;
d3_target      = sqrt(A_frontal_target*4/pi);
AR_diff_arr    = (d3_target/d2_init)^2;
AR_init        = AR_diff_arr;
theta_max_diff = 15;
flag           = "Friction";
d3             = d3_target;
m_dot_streamtube = Rho_1*v1*A2;
M_dot_in       = m_dot_streamtube;

[M3, P3, P3_0, A3, T3, T3_0, Rho_3, v3, M2, P2, P2_0, A2, T2, T2_0, Rho_2, v2, ...
    L_diffuser, M_dot_2, M_dot_3, d2, d3, Re_2, Re_3] = ...
    Diffuser_og_mine_for_7x7(M_dot_in, d3, M_dot_FOD, M_dot_comp, M1, P1, P1_0, ...
    T1, T1_0, Rho_1, v1, R, theta_max_diff, d2_init, A2, AR_init, a_inf, ...
    Gamma_inf, flag, mu_inf);

%%%--- SHARED: HX BOUNDARY CONDITIONS ---%%%

[T_mean_h, T_mean_c, T_c_i, T_c_o, T_h_i, T_h_o, C_h, C_c, C_star] = ...
    HX_deltaT(T3, M_dot_3, M_dot_coolant, T_out_fc, T_in_fc);

Pr_coolant = mu_EG(T_mean_h)*cp_EG_50_50(T_mean_h)/k_EG(T_mean_h);
Pr_air     = mu_air(T_mean_c)*cp_air(T_mean_c)/k_air(T_mean_c);


[dp_coolant,d_h_air, M_dot_4, b_t_air, b_t_coolant, dp_hx, N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, NTU, R_tot, v_channel_air, v_channel_coolant,d_h_coolant, A_o_coolant, A_o_air, Re_air, Re_coolant, h_air, h_coolant, L_solution, UA_unit, v4, P4_0, M4, T4, T4_0, P4, F_drag_hx, M_hx, A4, drag_HX, Q_pred_sol, T_cool_out, Pr_air, f_air, Nu_air, inlet_dp] = HX_design1(e,r,hx_theta, counter, x ,A3,v3, R,P3,P3_0, h, d3, T_h_o, n_modules, T_h_i, T_c_i, T_c_o, T_mean_h, T_mean_c, Q_tot, C_h, C_c, C_star,Pr_air, Pr_coolant, M_dot_coolant, M_dot_3, T3);
coolant_pres_drop = dp_coolant;

%% ----------------------------------------------- Puller fan ---------------------------------------------- %%
M_dot_des = rho_air(P4,T4)*A4*M4*sqrt(gamma_air(cp_air(T4))*R*T4); %design point mass flow rate
FPR_des = 1.01; %design point FPR
if fan == "ON"
    fpr_init = fan_fpr_simple(M_dot_4, M_dot_des, FPR_des);
    fprintf("Fan pressure ratio = %f\n",fpr_init);
else
    fpr_init = 1;
end

[M5,T5_0,T5,P5,P5_0,v5,A5,dp_fan,P_shaft,M_dot_5] = fan_backup (fpr_init,d3,P4,P4_0,T4,T4_0,v4,M4,A4,R,M_dot_4);

%% ----------------------------------------------- Nozzle design ----------------------------------------------%%%

% AR_noz = 0.33; % Nozzle area ratio
[A6, M6, L_nozzle, P6, P6_0, T6, T6_0, v6, M_dot_6, Re_5,d6,mdot_choked] = nozzle_old (AR_noz, flag,d3,P5,P5_0,T5,T5_0,v5,M5,M_dot_4);


%% ----------------------------------------------- Pressure difference correction ----------------------------------------------%%

% Without HX_deltaT in iteration
[P1,P2,P3,P4,P5,P6,P1_0,P2_0,P3_0,P4_0,P5_0,P6_0, M1, M2, M3, M4, M5, M6,  v1,v2,v3,v4,v5,v6, T1, T1_0, T2, T2_0, T3,T3_0,T4, T4_0, T5, T5_0, T6, T6_0, M_dot_2, M_dot_3, M_dot_4, M_dot_5, M_dot_6,P_shaft, M_hx, dp_coolant, L_solution, Q_pred_sol, T_cool_out, dp_hx, drag_HX, h_air, v_channel_air, Re_air, Pr_air, f_air, Nu_air, f_hx, inlet_dp] = pressure_eq_mDot_adjust_9thJan_noExhaust(e,r,m_dot_streamtube, fan, M_dot_des, FPR_des, hx_theta,fpr_init, M_dot_in, flight_phase, theta_max_diff,d3,M_dot_FOD, M_dot_comp,counter, A2_init, d2_init,AR_noz, dp_tot_prop, P1,P_inf, P_inf_tot, T_inf, V_inf, Rho_inf, AR_init, a_inf, Gamma_inf, flag, mu_inf,A2,v2, M2, R,P2,P2_0,h,d2,T2_0, T_h_o, n_modules, T_h_i, T_c_i, T_c_o, T_mean_h, T_mean_c, Q_tot, C_h, C_c, C_star,Pr_air, Pr_coolant, M_dot_coolant, M_dot_2, M_dot_3, T2,P3,P3_0,T3,T3_0,v3,M3,A4,P4,P4_0, T4, T4_0, v4, M4,P5,P5_0,M5,T5,T5_0,v5,P6,P6_0,M6,T6,T6_0,v6, M1, T1, T1_0, Rho_1, v1, P1_0, dia_prop);
% % With HX_deltaT in iteration
% [P1,P2,P3,P4,P5,P6,P1_0,P2_0,P3

if fpr_init == 1
    fprintf("The puller fan is not present\n")
else
    fprintf("The puller fan is present\n")
end


m_dot_spill = m_dot_streamtube-M_dot_2;
fprintf("Spillage mass flow rate = %f\n", m_dot_spill);
drag_tot = (m_dot_streamtube*v1) - (M_dot_6*v6) - (M_dot_2*M_dot_FOD*v3) + ((P2-P_inf)*A2) - ((P6-P_inf)*A6);


fprintf("Total drag = %f\n",drag_tot);

t_compute = toc(t_compute_start);
fprintf('Compute time: %.3f s\n', t_compute);

results.L_solution   = L_solution;
results.dp_hx        = dp_hx;
results.dp_coolant   = dp_coolant;
results.T4           = T4;
results.M_hx         = M_hx;
results.drag_HX      = drag_HX;
results.drag_tot     = drag_tot;
results.T_cool_out   = T_cool_out;
results.P_inf        = P_inf;
results.P1           = P1;
results.P2           = P2;
results.P3           = P3;
results.P4           = P4;
results.P5           = P5;
results.P6           = P6;
results.T_inf        = T_inf;
results.T1           = T1;
results.T2           = T2;
results.T3           = T3;
results.T5           = T5;
results.T6           = T6;
results.P_inf_tot    = P_inf_tot;
results.P1_0         = P1_0;
results.P2_0         = P2_0;
results.P3_0         = P3_0;
results.P4_0         = P4_0;
results.P5_0         = P5_0;
results.P6_0         = P6_0;
results.V_inf        = V_inf;
results.v1           = v1;
results.v2           = v2;
results.v3           = v3;
results.v4           = v4;
results.v5           = v5;
results.v6           = v6;
results.M_inf        = M_inf;
results.M1           = M1;
results.M2           = M2;
results.M3           = M3;
results.M4           = M4;
results.M5           = M5;
results.M6           = M6;
results.M_dot_2      = M_dot_2;
results.M_dot_3      = M_dot_3;
results.M_dot_4      = M_dot_4;
results.M_dot_5      = M_dot_5;
results.M_dot_6      = M_dot_6;
results.m_spill      = m_dot_streamtube - M_dot_2;
results.AR_init      = AR_init;
results.AR_noz       = AR_noz;
results.fpr          = fpr_init;
results.T_c_i        = T_c_i;
results.T_h_i        = T_h_i;
results.T_h_o        = T_h_o;
results.N_segments   = 1;
results.solve_T      = solve_T;
results.t_compute    = t_compute;
results.Q_pred_solution = Q_pred_sol;
results.h_air        = h_air;
results.v_channel_air= v_channel_air;
results.f_air        = f_air;
results.Re_air       = Re_air;
results.Pr_air       = Pr_air;
results.Nu_air       = Nu_air;
results.T_mean_c     = T_mean_c;
results.f_hx         = f_hx;
results.inlet_dp     = inlet_dp;
results.T_mean_h     = T_mean_h;

%%%--- SAVE (same scheme as run_disc_model_fwdpass; "lumped" in place of
%%%     the <N_segments>[-<N_cool_seg>][_dns] model-variant tag) ---%%%
if save_results

    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end

    % Physical input configuration -- built by the shared helper so a
    % caller doing its own cache lookup (via find_matching_id) before
    % calling this function builds the exact same key. Note fpr_init_arg
    % (captured before the fan-adjustment reassignment above), not
    % fpr_init, is passed here.
    inputs = physical_inputs(e, r, hx_theta, fan, fpr_init_arg, ...
        M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
        d2_init, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, ...
        solve_T);

    if exist(index_file, 'file')
        loaded = load(index_file, 'results_index');
        results_index = loaded.results_index;
    else
        results_index = table('Size', [0 2], 'VariableTypes', {'double', 'cell'}, ...
            'VariableNames', {'id', 'inputs'});
    end

    % Match against existing entries so the same physical-input case
    % always reuses the same id, whether it was first saved by this
    % function or by run_disc_model_fwdpass.
    id = [];
    for row = 1:height(results_index)
        if isequal(results_index.inputs{row}, inputs)
            id = results_index.id(row);
            break
        end
    end

    if isempty(id)
        id = height(results_index) + 1;
        new_row = table(id, {inputs}, 'VariableNames', {'id', 'inputs'});
        results_index = [results_index; new_row]; %#ok<AGROW>
        save(index_file, 'results_index');
    end

    base_fname = lumped_results_filename(id, solve_T);

    if overwrite
        results_file = fullfile(results_dir, base_fname);
    else
        results_file = fullfile(results_dir, ...
            next_versioned_filename(results_dir, base_fname));
    end

    save(results_file, 'results');
    fprintf('Saved results: %s\n', results_file);
end

end

% next_versioned_filename now lives in its own file (shared with
% run_disc_model_fwdpass.m), rather than as a local function here -- see
% next_versioned_filename.m.