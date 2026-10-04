clear all
clc

%%%--- SHARED INPUTS ---%%%

flight_phase = 3;

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

p11           = 22632;
t11           = 216.65;
R             = 287;
hx_theta      = 60;
fan           = "OFF";
e             = 4 - 2*0.1/sqrt(2); %6; %9.656; %10;
r             = 4 - 2*0.1/sqrt(2); %4.2; %3.5213; %4.2;
% n_modules     = 2;
d2_init       = 0.48; %0.52; % 0.46999; %0.44;
AR_diff       = 3.8; %3.6022; %4;
AR_noz        = 0.36; %0.38346; %0.33;
fpr_init      = 1;
M_dot_coolant = 44.4; % Max: 44.4, Min: 29
N_segments    = 50;

M_dot_FOD  = 0;
M_dot_comp = 0;

% d2_init = 0.45; %m, diffuser inlet initial guess
A2 = pi*d2_init*d2_init/4;
A_frontal_target = A2*AR_diff; %this number multiplied to A2 is your diffuser AR
d3_target = sqrt(A_frontal_target.*4/pi);
AR_diff_arr = (d3_target./d2_init).^2;
n_modules = 2; %ducts per nacelle

% tol_T = results_L.T_cool_out - results_L.T_h_o;
tol_T = 3;
use_DNS = false;
N_cool_seg = 0;   % old iterative 1D for this comparison; set >0 to compare against the 2D path instead

solve_T      = 0;
save_results = true;

N_list = [N_segments];

figs = struct();
col_D   = autumn(length(N_list));
col_DNS = winter(length(N_list));

%%%--- LUMPED MODEL -- check cache, run only if needed ---%%%
% Same pattern as main_disc_test.m: build the matching key with
% physical_inputs, look it up (read-only) with find_matching_id, and only
% load instead of running if BOTH an id match and the actual results file
% for this variant exist.

inputs_L = physical_inputs(e, r, hx_theta, fan, fpr_init, M_dot_coolant, n_modules, ...
    Q_tot, T_in_fc, T_out_fc, h, p11, t11, d2_init, AR_diff, AR_noz, V_inf, R, ...
    flight_phase, M_dot_FOD, M_dot_comp, solve_T);

id_L = find_matching_id(inputs_L);

results_L = [];
if ~isempty(id_L)
    candidate_file = fullfile('results', lumped_results_filename(id_L, solve_T));
    if exist(candidate_file, 'file')
        fprintf('Found cached results, loading: %s\n', candidate_file);
        loaded    = load(candidate_file, 'results');
        results_L = loaded.results;
    end
end

if isempty(results_L)
    results_L = run_lumped_model(e, r, hx_theta, fan, fpr_init, M_dot_coolant, n_modules, ...
        Q_tot, T_in_fc, T_out_fc, h, p11, t11, d2_init, AR_diff, AR_noz, V_inf, R, ...
        flight_phase, M_dot_FOD, M_dot_comp, solve_T, save_results);
end

figs = plot_stations(results_L, 'lumped', 'figs', figs);
figs = plot_HX(results_L, 'lumped', 'figs', figs);

for i = 1:length(N_list)

    %%%--- DISCRETISED MODEL -- check cache, run only if needed ---%%%

    inputs_D = physical_inputs(e, r, hx_theta, fan, fpr_init, M_dot_coolant, n_modules, ...
        Q_tot, T_in_fc, T_out_fc, h, p11, t11, d2_init, AR_diff, AR_noz, V_inf, R, ...
        flight_phase, M_dot_FOD, M_dot_comp, solve_T, tol_T);

    id_D = find_matching_id(inputs_D);

    results_D = [];
    if ~isempty(id_D)
        candidate_file = fullfile('results', ...
            disc_results_filename(N_list(i), N_cool_seg, false, id_D, solve_T));
        if exist(candidate_file, 'file')
            fprintf('Found cached results, loading: %s\n', candidate_file);
            loaded    = load(candidate_file, 'results');
            results_D = loaded.results;
        end
    end

    if isempty(results_D)
        results_D = run_disc_model_fwdpass(false, N_list(i), N_cool_seg, e, r, hx_theta, fan, fpr_init, ...
            M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
            d2_init, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, ...
            solve_T, save_results, tol_T);
    end

    figs = plot_stations(results_D, 'discretised', 'figs', figs, 'color', col_D(i,:));
    figs = plot_HX(results_D, 'discretised', 'figs', figs, 'color', col_D(i,:));

    if use_DNS
        %%%--- DNS VARIANT -- same physical inputs/id as above, just a
        %%%     different results file (use_DNS = true) ---%%%
        results_DNS = [];
        if ~isempty(id_D)
            candidate_file_dns = fullfile('results', ...
                disc_results_filename(N_list(i), N_cool_seg, true, id_D, solve_T));
            if exist(candidate_file_dns, 'file')
                fprintf('Found cached results, loading: %s\n', candidate_file_dns);
                loaded      = load(candidate_file_dns, 'results');
                results_DNS = loaded.results;
            end
        end

        if isempty(results_DNS)
            results_DNS = run_disc_model_fwdpass(true, N_list(i), N_cool_seg, e, r, hx_theta, fan, fpr_init, ...
                M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
                d2_init, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, ...
                solve_T, save_results, tol_T);
        end

        figs = plot_stations(results_DNS, 'dns', 'figs', figs, 'color', col_DNS(i,:));
        figs = plot_HX(results_DNS, 'dns', 'figs', figs, 'color', col_DNS(i,:));
    end
end

results_cfd = process_wall_CFD_data(...
    "C:\Users\tiesh\TUDelft\Thesis\CFD\\4by4_150ch_wall_q_air", ...
    "C:\Users\tiesh\TUDelft\Thesis\CFD\\4by4_150ch_wall_Twall_air", ...
    "C:\Users\tiesh\TUDelft\Thesis\CFD\\4by4_150ch_wall_Tbulk_air", ...
    results_D.d_h_air);
% figs = plot_HX(results_cfd, 'cfd',         'figs', figs);



% fprintf("\n========================================\n");
% fprintf("           MODEL COMPARISON             \n");
% fprintf("========================================\n");
% fprintf("%-30s %-15s %-15s\n", "Quantity", "Lumped", "Discretised");
% fprintf("%-30s %-15.4f %-15.4f\n", "HX length [m]",       results_L.L_solution,  results_D.L_solution);
% fprintf("%-30s %-15.2f %-15.2f\n", "HX dp [Pa]",          results_L.dp_hx,       results_D.dp_hx);
% fprintf("%-30s %-15.2f %-15.2f\n", "Air outlet T4 [K]",   results_L.T4,          results_D.T4);
% fprintf("%-30s %-15.2f %-15.2f\n", "HX mass [kg]",        results_L.M_hx,        results_D.M_hx);
% fprintf("%-30s %-15.2f %-15.2f\n", "HX drag [N]",         results_L.drag_HX,     results_D.drag_HX);
% fprintf("%-30s %-15.2f %-15.2f\n", "P6 - P_inf [Pa]",     results_L.P6-results_L.P_inf,    results_D.P6-results_D.P_inf);
% fprintf("%-30s %-15.2f %-15.2f\n", "Coolant outlet temp [K]",  results_L.T_cool_out,  results_D.T_cool_out);
% fprintf("%-30s %-15.2f %-15.2f\n", "Coolant pressure drop [bar]",  results_L.dp_coolant/1e5,  results_D.dp_coolant/1e5);
% fprintf("%-30s %-15.4f %-15.4f\n", "Cooling power [MW]", results_L.Q_pred_solution/1e6,     results_D.Q_pred_solution/1e6);
% fprintf("%-30s %-15.4f %-15.4f\n", "Mass flow in [kg/s]", results_L.M_dot_2,     results_D.M_dot_2);
% fprintf("%-30s %-15.4f %-15.4f\n", "Spillage [kg/s]",     results_L.m_spill,     results_D.m_spill);
% fprintf("%-30s %-15.2f %-15.2f\n", "Bulk inlet pressure drag [Pa]",     results_L.inlet_dp,     results_D.inlet_dp);