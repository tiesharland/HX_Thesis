function results = run_disc_model_fwdpass(use_DNS, N_segments, N_cool_seg, e, r, hx_theta, fan, fpr_init, ...
    M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
    d2_init, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, tol_T)

counter = 1;

%%%--- FREESTREAM ---%%%
[P_inf, T_inf, Rho_inf, Gamma_inf, Cp_inf, mu_inf, P_inf_tot, M_inf, a_inf] = ...
    Freestream(h, p11, R, t11, V_inf);
P_inf_tot = P_inf + 0.5*Rho_inf*V_inf^2;
a_inf     = sqrt(Gamma_inf*R*T_inf);
M_inf     = V_inf/a_inf;

%%%--- PROPELLER ---%%%
[J_val, Cp_val, eff_val, D, J, Cp, V_ind_prop, dp_tot_prop, T_net, dia_prop, ...
    M1, P1, P1_0, T1, T1_0, Rho_1, v1, Re_1] = ...
    propeller(flight_phase, V_inf, Rho_inf, P_inf, a_inf, T_inf, mu_inf);

%%%--- DIFFUSER GEOMETRY (fixed, computed once) ---%%%
A2               = pi*d2_init^2/4;
A_frontal_target = A2*AR_diff;
d3_target        = sqrt(A_frontal_target*4/pi);
AR_diff_arr      = (d3_target/d2_init)^2;
AR_init          = AR_diff_arr;
theta_max_diff   = 15;
flag             = "Friction";
d3               = d3_target;
m_dot_streamtube = Rho_1*v1*A2;

%%%--- SINGLE FORWARD PASS FUNCTION ---%%%
% Takes M_dot_in, returns full state and diff_P = P6 - P_inf.
% All upstream quantities (propeller, freestream, geometry) captured
% from the enclosing workspace via closure.

    function [diff_P, state] = forward_pass(M_dot_in)

        fprintf('M_dot_in = %.4f, m_dot_streamtube = %.4f\n', M_dot_in, m_dot_streamtube)

        %--- Diffuser ---%
        [M3, P3, P3_0, A3, T3, T3_0, Rho_3, v3, M2, P2, P2_0, A2_loc, T2, T2_0, ...
         Rho_2, v2, L_diffuser, M_dot_2, M_dot_3, d2, d3_loc, Re_2, Re_3] = ...
            Diffuser_og_mine_for_7x7(M_dot_in, d3, M_dot_FOD, M_dot_comp, M1, ...
            P1, P1_0, T1, T1_0, Rho_1, v1, R, theta_max_diff, d2_init, A2, ...
            AR_init, a_inf, Gamma_inf, flag, mu_inf);

        %--- HX boundary conditions ---%
        [T_mean_h, T_mean_c, T_c_i, T_c_o, T_h_i, T_h_o, C_h, C_c, C_star] = ...
            HX_deltaT(T3, M_dot_3, M_dot_coolant, T_out_fc, T_in_fc);

        %--- Discretised HX ---%

        [dp_coolant_loop, d_h_air, M_dot_4, b_t_air, b_t_coolant, dp_hx, ...
            N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, NTU, R_tot, ...
            v_channel_air, v_channel_coolant, d_h_coolant, A_o_coolant, A_o_air, ...
            Re_air, Re_coolant, h_air, h_coolant, L_solution, UA_unit, v4, P4_0, ...
            M4, T4, T4_0, P4, F_drag_hx, M_hx, A4, drag_HX, T_cool_seg, dp_cool_seg, ...
            T_air_seg, P_air_seg, v_air_seg, Re_air_seg, Pr_air_seg, v_channel_seg, ...
            K_seg, f_air_seg, Nu_air_seg, h_air_seg, eta_fin_seg, ...
            UA_seg_arr, NTU_seg_arr, eps_seg_arr, Q_seg_arr, dp_seg_arr, ...
            Q_pred_solution, T_h_o_solution, P_0_air_seg, T_0_air_seg, M_air_seg, f_hx_seg, inlet_dp, outlet_dp, ...
            delta_BL_seg, A_free_seg, d_h_bulk_seg, T_mean_c_arr, T_mean_h_arr] = ...
            HX_design1_disc(use_DNS, e, r, hx_theta, counter, A3, v3, R, P3, ...
            d3_loc, T_h_o, n_modules, T_h_i, T_c_i, T_mean_h, ...
            Q_tot, M_dot_coolant, M_dot_3, T3, N_segments, N_cool_seg, tol_T);

        % if N_cool_seg
        % 
        %     [dp_coolant_loop, d_h_air, M_dot_4, b_t_air, b_t_coolant, dp_hx, ...
        %         N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, NTU, R_tot, ...
        %         v_channel_air, v_channel_coolant, d_h_coolant, A_o_coolant, A_o_air, ...
        %         Re_air, Re_coolant, h_air, h_coolant, L_solution, UA_unit, v4, P4_0, ...
        %         M4, T4, T4_0, P4, F_drag_hx, M_hx, A4, drag_HX, T_cool_seg, dp_cool_seg, ...
        %         T_air_seg, P_air_seg, v_air_seg, Re_air_seg, Pr_air_seg, v_channel_seg, ...
        %         K_seg, f_air_seg, Nu_air_seg, h_air_seg, eta_fin_seg, ...
        %         UA_seg_arr, NTU_seg_arr, eps_seg_arr, Q_seg_arr, dp_seg_arr, ...
        %         Q_pred_solution, T_h_o_solution, P_0_air_seg, T_0_air_seg, M_air_seg, f_hx_seg, inlet_dp, outlet_dp, ...
        %         delta_BL_seg, A_free_seg, d_h_bulk_seg, T_mean_c_arr, T_mean_h_arr] = ...
        %         HX_design1_disc_2d_v2(use_DNS, e, r, hx_theta, counter, A3, v3, R, P3, ...
        %         d3_loc, T_h_o, n_modules, T_h_i, T_c_i, T_mean_h, ...
        %         Q_tot, M_dot_coolant, M_dot_3, T3, N_segments, N_cool_seg, tol_T);
        % else
        % 
        %     [dp_coolant_loop, d_h_air, M_dot_4, b_t_air, b_t_coolant, dp_hx, ...
        %         N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, NTU, R_tot, ...
        %         v_channel_air, v_channel_coolant, d_h_coolant, A_o_coolant, A_o_air, ...
        %         Re_air, Re_coolant, h_air, h_coolant, L_solution, UA_unit, v4, P4_0, ...
        %         M4, T4, T4_0, P4, F_drag_hx, M_hx, A4, drag_HX, T_cool_seg, dp_cool_seg, ...
        %         T_air_seg, P_air_seg, v_air_seg, Re_air_seg, Pr_air_seg, v_channel_seg, ...
        %         K_seg, f_air_seg, Nu_air_seg, h_air_seg, eta_fin_seg, ...
        %         UA_seg_arr, NTU_seg_arr, eps_seg_arr, Q_seg_arr, dp_seg_arr, ...
        %         Q_pred_solution, T_h_o_solution, P_0_air_seg, T_0_air_seg, M_air_seg, f_hx_seg, inlet_dp, outlet_dp, ...
        %         delta_BL_seg, A_free_seg, d_h_bulk_seg, T_mean_c_arr, T_mean_h_arr] = ...
        %         HX_design1_disc(use_DNS, e, r, hx_theta, counter, A3, v3, R, P3, ...
        %         d3_loc, T_h_o, n_modules, T_h_i, T_c_i, T_mean_h, ...
        %         Q_tot, M_dot_coolant, M_dot_3, T3, N_segments, tol_T);
        % 
        % end

        %--- Fan ---%
        fpr_loc = 1;
        [M5, T5_0, T5, P5, P5_0, v5, A5, dp_fan, P_shaft, M_dot_5] = ...
            fan_backup(fpr_loc, d3_loc, P4, P4_0, T4, T4_0, v4, M4, A4, R, M_dot_4);

        %--- Nozzle ---%
        [A6, M6, L_nozzle, P6, P6_0, T6, T6_0, v6, M_dot_6, ~, d6, ~] = ...
            nozzle_old(AR_noz, flag, d3_loc, P5, P5_0, T5, T5_0, v5, M5, M_dot_4);

        diff_P = P6 - P_inf;
        fprintf("diff_P = %.4f  M_dot_in = %.4f\n", diff_P, M_dot_in);

        %--- Pack state ---%
        state.P1 = P1; state.P2 = P2; state.P3 = P3;
        state.P4 = P4; state.P5 = P5; state.P6 = P6;
        state.P1_0 = P1_0; state.P2_0 = P2_0; state.P3_0 = P3_0;
        state.P4_0 = P4_0; state.P5_0 = P5_0; state.P6_0 = P6_0;
        state.T1 = T1; state.T2 = T2; state.T3 = T3;
        state.T4 = T4; state.T5 = T5; state.T6 = T6;
        state.T1_0 = T1_0; state.T2_0 = T2_0; state.T3_0 = T3_0;
        state.T4_0 = T4_0; state.T5_0 = T5_0; state.T6_0 = T6_0;
        state.M1 = M1; state.M2 = M2; state.M3 = M3;
        state.M4 = M4; state.M5 = M5; state.M6 = M6;
        state.v1 = v1; state.v2 = v2; state.v3 = v3;
        state.v4 = v4; state.v5 = v5; state.v6 = v6;
        state.M_dot_2 = M_dot_2; state.M_dot_3 = M_dot_3;
        state.M_dot_4 = M_dot_4; state.M_dot_5 = M_dot_5;
        state.M_dot_6 = M_dot_6;
        state.P_shaft = P_shaft;
        state.M_hx = M_hx; state.drag_HX = drag_HX;
        state.dp_hx = dp_hx; state.dp_coolant = dp_coolant_loop;
        state.L_solution = L_solution;
        state.Q_pred_solution = Q_pred_solution;
        state.T_h_o_solution = T_h_o_solution;
        state.T_air_seg = T_air_seg; state.P_air_seg = P_air_seg;
        state.v_air_seg = v_air_seg; state.v_channel_seg = v_channel_seg;
        state.Re_air_seg = Re_air_seg; state.Pr_air_seg = Pr_air_seg;
        state.f_air_seg = f_air_seg; state.Nu_air_seg = Nu_air_seg;
        state.h_air_seg = h_air_seg;
        state.Q_seg_arr = Q_seg_arr; state.dp_seg_arr = dp_seg_arr;
        state.T_cool_seg = T_cool_seg;
        state.dp_cool_seg = dp_cool_seg;
        state.P_0_air_seg = P_0_air_seg; state.T_0_air_seg = T_0_air_seg;
        state.M_air_seg = M_air_seg;
        state.f_hx_seg = f_hx_seg; 
        state.inlet_dp = inlet_dp;
        state.outlet_dp = outlet_dp;
        state.A_o_air = A_o_air; 
        state.N_fin_air = N_fin_air; state.N_fin_coolant = N_fin_coolant;
        state.N_air_pass = N_air_pass; state.N_coolant_pass = N_coolant_pass;
        state.T_h_i = T_h_i; state.T_h_o = T_h_o; state.T_c_i = T_c_i;
        state.d_h_air = d_h_air;
        state.delta_BL_seg = delta_BL_seg;
        state.A_free_seg = A_free_seg;
        state.d_h_bulk_seg = d_h_bulk_seg;
        state.T_mean_c_arr = T_mean_c_arr; state.T_mean_h_arr = T_mean_h_arr;

    end

%%%--- PRESSURE BALANCE ITERATION ---%%%

% p_loop   = P_inf;
% diff_P   = Inf;
% started  = false;

M_dot_in = m_dot_streamtube;

% Turning point detection and initialisation

[diff_P, state] = forward_pass(M_dot_in);
fprintf("diff_P = %.4f  M_dot_in = %.4f\n", diff_P, M_dot_in);

diff_P_prev   = diff_P;   % track previous diff_P to detect worsening
diff_P_best   = diff_P;  % track least-negative diff_P seen so far
M_dot_in_best = M_dot_in;
state_best      = state;

while abs(diff_P) > 10% || ~started
    % started = true;

    % Detect turning point: diff_P worsening after previously improving
    if  abs(diff_P) > abs(diff_P_prev)
        fprintf(['Warning: diff_P worsened (%.4f -> %.4f). ' ...
                 'Turning point passed.\n'], diff_P_prev, diff_P);
        fprintf('Best achievable: diff_P = %.4f at M_dot_in = %.4f\n', ...
                diff_P_best, M_dot_in_best);

        if abs(diff_P_best) < 10
            fprintf('Within 10 Pa tolerance — accepting best point.\n');
            diff_P = diff_P_best;
            state    = state_best;
            break
        else
            error(['No feasible pressure balance found. ' ...
                   'Best diff_P = %.4f Pa at M_dot_in = %.4f kg/s.\n' ...
                   'Design infeasible — consider increasing d2_init, ' ...
                   'reducing e or r, or increasing M_dot_coolant.'], ...
                   diff_P_best, M_dot_in_best);
        end
    end

    diff_P_prev = diff_P;

    %--- Same step-based mass flow correction as pressure_iteration_disc ---%
    if diff_P < 0
        if diff_P < -3000
            M_dot_in = M_dot_in - (0.1*state.M_dot_3);
        elseif diff_P > -3000 && diff_P < -1000
            M_dot_in = M_dot_in - (0.01*state.M_dot_3);
        elseif diff_P > -1000 && diff_P < -200
            M_dot_in = M_dot_in - (0.005*state.M_dot_3);
        else
            M_dot_in = M_dot_in - (0.0005*state.M_dot_3);
        end
    elseif diff_P > 0
        if diff_P > 1000
            M_dot_in = M_dot_in + (0.01*state.M_dot_3);
        elseif diff_P < 1000 && diff_P > 200
            M_dot_in = M_dot_in + (0.005*state.M_dot_3);
        else
            M_dot_in = M_dot_in + (0.0005*state.M_dot_3);
        end
    else
        break
    end

    if fan == "OFF" && M_dot_in > m_dot_streamtube
        error("Mass flow rate into diffuser is not sufficient; Puller fan is needed;")
    end
    
    [diff_P, state] = forward_pass(M_dot_in);
    fprintf("diff_P = %.4f  M_dot_in = %.4f\n", diff_P, M_dot_in);

    % fprintf('P1 = %.2f, P6 = %.2f, diff_P = %.6e, m_dot = %.4f\n', ...
    %         P1, state.P6, diff_P, M_dot_in);
    
    % Track best (least negative) operating point
    if abs(diff_P) < abs(diff_P_best)
        diff_P_best   = diff_P;
        M_dot_in_best = M_dot_in;
        state_best      = state;
    end

end

fprintf('Converged: |P6 - P_inf| = %.4f Pa at M_dot_in = %.4f kg/s\n', ...
        abs(diff_P), M_dot_in);

%%%--- COMPUTE m_dot_seg ---%%%
m_dot_seg = zeros(1, N_segments+1);
for k = 1:N_segments+1
    rho_k        = rho_air(state.P_air_seg(k), state.T_air_seg(k));
    m_dot_seg(k) = rho_k * state.v_channel_seg(k) * state.A_o_air * ...
                   state.N_fin_air * state.N_air_pass;
end
m_dot_seg(1) = m_dot_seg(2);

%%%--- PACK RESULTS ---%%%
m_dot_spill = m_dot_streamtube - state.M_dot_2;
A6 = AR_noz * pi*d3^2/4;
drag_tot = (m_dot_streamtube*v1) - (state.M_dot_6*state.v6) ...
         - (state.M_dot_2*M_dot_FOD*state.v3) ...
         + (state.P2-P_inf)*A2 - (state.P6-P_inf)*A6;


results.L_solution      = state.L_solution;
results.dp_hx           = state.dp_hx;
results.dp_coolant      = state.dp_coolant;
results.T4              = state.T4;
results.M_hx            = state.M_hx;
results.drag_HX         = state.drag_HX;
results.drag_tot        = drag_tot;
results.T_cool_out      = state.T_h_o_solution;
results.P_inf           = P_inf;
results.P_inf_tot       = P_inf_tot;
results.P1              = state.P1;   results.P2 = state.P2;
results.P3              = state.P3;   results.P4 = state.P4;
results.P5              = state.P5;   results.P6 = state.P6;
results.P1_0            = state.P1_0; results.P2_0 = state.P2_0;
results.P3_0            = state.P3_0; results.P4_0 = state.P4_0;
results.P5_0            = state.P5_0; results.P6_0 = state.P6_0;
results.T_inf           = T_inf;
results.T1              = state.T1;   results.T2 = state.T2;
results.T3              = state.T3;   results.T5 = state.T5;
results.T6              = state.T6;
results.V_inf           = V_inf;
results.v1              = v1;
results.v2              = state.v2;   results.v3 = state.v3;
results.v4              = state.v4;   results.v5 = state.v5;
results.v6              = state.v6;
results.M_inf           = M_inf;
results.M1              = M1;
results.M2              = state.M2;   results.M3 = state.M3;
results.M4              = state.M4;   results.M5 = state.M5;
results.M6              = state.M6;
results.M_dot_2         = state.M_dot_2;
results.M_dot_3         = state.M_dot_3;
results.M_dot_4         = state.M_dot_4;
results.M_dot_5         = state.M_dot_5;
results.M_dot_6         = state.M_dot_6;
results.m_spill         = m_dot_spill;
results.m_dot_seg       = m_dot_seg;
results.AR_init         = AR_init;
results.AR_noz          = AR_noz;
results.fpr             = fpr_init;
results.T_c_i           = state.T_c_i;
results.T_h_i           = state.T_h_i;
results.T_h_o           = state.T_h_o;
results.N_segments      = N_segments;
results.T_air_seg       = state.T_air_seg;
results.P_air_seg       = state.P_air_seg;
results.v_air_seg       = state.v_air_seg;
results.v_channel_seg   = state.v_channel_seg;
results.Re_air_seg      = state.Re_air_seg;
results.Pr_air_seg      = state.Pr_air_seg;
results.f_air_seg       = state.f_air_seg;
results.Nu_air_seg      = state.Nu_air_seg;
results.h_air_seg       = state.h_air_seg;
results.Q_seg_arr       = state.Q_seg_arr;
results.dp_seg_arr      = state.dp_seg_arr;
results.T_cool_seg      = state.T_cool_seg;
results.dp_cool_seg     = state.dp_cool_seg;
results.Q_pred_solution = state.Q_pred_solution;
results.P_0_air_seg     = state.P_0_air_seg;
results.T_0_air_seg     = state.T_0_air_seg;
results.M_air_seg       = state.M_air_seg;
results.f_hx_seg        = state.f_hx_seg;
results.inlet_dp        = state.inlet_dp;
results.outlet_dp       = state.outlet_dp;
results.d_h_air         = state.d_h_air;
results.delta_BL_seg   = state.delta_BL_seg;
results.A_free_seg      = state.A_free_seg;
results.d_h_bulk_seg    = state.d_h_bulk_seg;
results.T_mean_h_arr = state.T_mean_h_arr; 
results.T_mean_c_arr = state.T_mean_c_arr;
results.A_o_air = state.A_o_air; 
results.N_fin_air = state.N_fin_air; results.N_fin_coolant = state.N_fin_coolant;
results.N_air_pass = state.N_air_pass; results.N_coolant_pass = state.N_coolant_pass;
end