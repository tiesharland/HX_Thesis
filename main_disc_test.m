clear all
clc

counter = 1;

e = 4 - 2*0.1/sqrt(2);
r = 4 - 2*0.1/sqrt(2);
d2_init = 0.48; %0.45808;
AR_diff = 3.8;
AR_noz  = 0.36;
A2      = pi*d2_init*d2_init/4;
A_frontal_target = A2*AR_diff;
d3_target        = sqrt(A_frontal_target*4/pi);
AR_diff_arr      = (d3_target/d2_init)^2;
n_modules        = 2;
M_dot_coolant    = 44.4;

%%%--- INPUTS - FREESTREAM ---%%%

flight_phase = 3; % 1=take-off, 2=top of climb, 3=cruise
if flight_phase == 1
    h        = 5;
    T_in_fc  = 85+273;
    T_out_fc = 105+273;
    Q_tot    = 2.38e6;
    V_inf    = 2;
elseif flight_phase == 3
    h        = 4876;
    T_in_fc  = 70+273;
    T_out_fc = 85+273;
    V_inf    = 128;
    Q_tot    = 2.25e6;
elseif flight_phase == 2
    h        = 4876;
    T_in_fc  = 70+273;
    T_out_fc = 85+273;
    Q_tot    = 2.25e6;
    V_inf    = 128;
end

p11    = 22632;
t11    = 216.65;
R      = 287;
hx_theta = 60;
fpr_init = 1;

%%%--- MASS FLOW LEAKS ---%%%

M_dot_FOD  = 0;
M_dot_comp = 0;

%%%--- DIFFUSER INPUTS ---%%%

AR_init = AR_diff_arr;

%%%--- RUN MODEL ---%%%

fan  = "OFF";
disc = "ON";
N_segments = 50;
N_cool_seg = 0;
tol_T = 0; % K
use_DNS = true;
solve_T = 0;
save_results = true;
overwrite = false;

% Set true to always re-run the model, bypassing the cache lookup below
% entirely (e.g. to get a fresh t_compute for comparison after a
% performance change). The run still saves per save_results/overwrite as
% usual, so a false->true->false round trip leaves the cache consistent.
force_recompute = true;

%%%--- CHECK CACHE, RUN ONLY IF NEEDED ---%%%
% 1. Build the same "physical input configuration" key
%    run_disc_model_fwdpass itself would build when saving.
% 2. Look it up in results/results_index.mat (find_matching_id does NOT
%    create or modify the index -- it's a read-only lookup).
% 3. An id match only means SOME variant at these physical inputs was
%    saved before -- it says nothing about whether THIS variant
%    (N_segments, N_cool_seg, use_DNS) was. So the actual file for this
%    variant+id still has to be checked before trusting it's there.
% 4. Only if both the id and that specific file exist do we load instead
%    of running.

inputs = physical_inputs(e, r, hx_theta, fan, fpr_init, ...
    M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
    d2_init, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, ...
    solve_T, tol_T);

id = find_matching_id(inputs);

results_D        = [];
loaded_from_cache = false;

if ~isempty(id) && ~force_recompute
    candidate_file = fullfile('results', ...
        disc_results_filename(N_segments, N_cool_seg, use_DNS, id, solve_T));
    if exist(candidate_file, 'file')
        fprintf('Found cached results, loading: %s\n', candidate_file);
        loaded    = load(candidate_file, 'results');
        results_D = loaded.results;
        loaded_from_cache = true;
    end
elseif ~isempty(id) && force_recompute
    fprintf('force_recompute is true -- skipping cache and re-running.\n');
end

if isempty(results_D)
    % Fresh DNS queries only come from an actual run, so the coverage log
    % only means something right before the model executes -- reset it
    % here, not earlier (a cache hit never reaches this branch).
    if use_DNS
        thermoturb_query_log('reset');
    end

    results_D = run_disc_model_fwdpass(use_DNS, N_segments, N_cool_seg, e, r, hx_theta, fan, fpr_init, ...
        M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
        d2_init, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, ...
        solve_T, save_results, overwrite, tol_T);
end

%%%--- SUMMARY PRINTS ---%%%
fprintf("\n=== HX RESULTS ===\n");
fprintf("HX length                      = %.4f m\n",   results_D.L_solution);
fprintf("HX pressure drop               = %.2f Pa\n",  results_D.dp_hx);
fprintf("Air outlet temp                = %.2f K\n",   results_D.T4);
fprintf("HX mass                        = %.2f kg\n",  results_D.M_hx);
fprintf("HX drag                        = %.2f N\n",   results_D.drag_HX);
fprintf("P6 - P_inf                     = %.2f Pa\n",  results_D.P6 - results_D.P_inf);
fprintf("Coolant outlet temp            = %.2f K\n",   results_D.T_cool_out);
fprintf("Coolant pressure drop (avg)    = %.2f Pa\n",   results_D.dp_coolant);
fprintf("Cooling power                  = %.3f MW\n",   results_D.Q_pred_solution/1e6);
fprintf("Nozzle exit P6                 = %.2f Pa\n",  results_D.P6);
fprintf("Ambient P_inf                  = %.2f Pa\n",  results_D.P_inf);
fprintf("d_h / L                        = %f \n", results_D.d_h_air/results_D.L_solution);

%%%--- PLOTS ---%%%

figs = plot_stations(results_D, 'discretised');
figs = plot_HX(results_D, 'discretised', 'figs', figs);

if use_DNS
    if loaded_from_cache
        fprintf(['Skipping DNS coverage plot: results were loaded from cache, ' ...
                 'so no DNS queries were logged this run.\n']);
    else
        plot_thermoturb_coverage();
    end
end

% results_cfd = process_wall_CFD_data(...
%     "C:\Users\tiesh\TUDelft\Thesis\CFD\4by4_150ch_wall\run 2\q_air", ...
%     "C:\Users\tiesh\TUDelft\Thesis\CFD\4by4_150ch_wall\run 2\Twall_air", ...
%     "C:\Users\tiesh\TUDelft\Thesis\CFD\4by4_150ch_wall\run 2\Tbulk_air", ...
%     "C:\Users\tiesh\TUDelft\Thesis\CFD\4by4_150ch_wall\run 2\Pbulk_air", ...
%     "C:\Users\tiesh\TUDelft\Thesis\CFD\4by4_150ch_wall\run 2\Vmag_air", ...
%     "C:\Users\tiesh\TUDelft\Thesis\CFD\4by4_150ch_wall\run 2\M_dot_air", ...
%     results_D.d_h_air);
%
% compare_cfd_vs_2d(results_D, results_cfd);