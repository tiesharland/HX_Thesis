%% compare_models.m
% Master script for comparing models on ONE input case with plot_stations
% and plot_HX. Set the case below, then choose what to compare:
%
%   include_lumped  - add the lumped model
%   use_DNS_list    - correlation-based (false) and/or DNS-based (true)
%   N_cool_list     - 0 = old iterative 1D, >= 1 = 2D direct marching with
%                     that many coolant-side segments
%   N_list          - air-side segments
%
% Every combination of use_DNS_list x N_cool_list x N_list is one line in
% the figures. Colour = (DNS or not, N_cool_seg); line style = N_segments.
% E.g. N_cool_list = [0 10], N_list = [10 50] gives four discretised lines.
%
% Each model evaluation first looks for a saved result for exactly these
% inputs (results/id-<id>/...) and only runs the model if none exists, as
% in main_disc_test.m. New runs are saved per save_results / overwrite.
% force_recompute = true skips the lookup and always runs (see below).

clear all
clc

%%%--- INPUT CASE ---%%%

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
d2_init       = 0.48; %0.52; % 0.46999; %0.44;
AR_diff       = 3.8; %3.6022; %4;
AR_noz        = 0.36; %0.38346; %0.33;
fpr_init      = 1;
M_dot_coolant = 44.4; % Max: 44.4, Min: 29
n_modules     = 2;    %ducts per nacelle

M_dot_FOD  = 0;
M_dot_comp = 0;

% tol_T = results_L.T_cool_out - results_L.T_h_o;
tol_T   = 0;
solve_T = 0;

%%%--- WHAT TO COMPARE ---%%%

include_lumped = true;
use_DNS_list   = false;   % [false true] to compare correlation-based vs DNS-based
N_cool_list    = 0;       % 0 = 1D iterative; e.g. [0 5 10] adds 2D lines
N_list         = 50;      % e.g. [10 50] to compare air-side resolutions

line_styles = {'-', '--', ':', '-.'};   % one per entry of N_list, in order

%%%--- RUN / CACHE OPTIONS ---%%%

save_results    = true;
overwrite       = false;   % keep existing cached files; save new versions instead of clobbering
force_recompute = false;   % true = ignore saved results and re-run everything

include_cfd = false;       % true = also process the CFD data (needs the CFD folders on this machine)

%% ---- Derived (kept from the old script) ---- %%
A2 = pi*d2_init*d2_init/4;
A_frontal_target = A2*AR_diff; %this number multiplied to A2 is your diffuser AR
d3_target = sqrt(A_frontal_target.*4/pi);
AR_diff_arr = (d3_target./d2_init).^2; %#ok<NASGU>

figs = struct();
summary = {};   % rows: label, L, dp_hx, T4, M_hx, drag_HX, from_cache

any_dns = any(use_DNS_list);
if any_dns
    thermoturb_query_log('reset');   % so the coverage plot (below) only shows this script's fresh runs
end
ran_dns_fresh = false;

%%%--- LUMPED MODEL ---%%%
% Lumped results are keyed without tol_T (see physical_inputs.m).

if include_lumped
    inputs_L = physical_inputs(e, r, hx_theta, fan, fpr_init, M_dot_coolant, n_modules, ...
        Q_tot, T_in_fc, T_out_fc, h, p11, t11, d2_init, AR_diff, AR_noz, V_inf, R, ...
        flight_phase, M_dot_FOD, M_dot_comp, solve_T);

    id_L = find_matching_id(inputs_L);

    results_L = [];
    from_cache = false;
    if ~isempty(id_L) && ~force_recompute
        candidate_file = fullfile('results', lumped_results_filename(id_L, solve_T));
        if exist(candidate_file, 'file')
            fprintf('Found cached results, loading: %s\n', candidate_file);
            loaded     = load(candidate_file, 'results');
            results_L  = loaded.results;
            from_cache = true;
        end
    end

    if isempty(results_L)
        results_L = run_lumped_model(e, r, hx_theta, fan, fpr_init, M_dot_coolant, n_modules, ...
            Q_tot, T_in_fc, T_out_fc, h, p11, t11, d2_init, AR_diff, AR_noz, V_inf, R, ...
            flight_phase, M_dot_FOD, M_dot_comp, solve_T, save_results);
    end

    figs = plot_stations(results_L, 'lumped', 'figs', figs);
    figs = plot_HX(results_L,       'lumped', 'figs', figs);
    summary(end+1,:) = {'Lumped', results_L.L_solution, results_L.dp_hx, results_L.T4, ...
                        results_L.M_hx, results_L.drag_HX, from_cache};
end

%%%--- DISCRETISED VARIANTS ---%%%
% Same physical inputs -> same id for every variant; only the results file
% differs (N_segments, N_cool_seg, use_DNS).

inputs_D = physical_inputs(e, r, hx_theta, fan, fpr_init, M_dot_coolant, n_modules, ...
    Q_tot, T_in_fc, T_out_fc, h, p11, t11, d2_init, AR_diff, AR_noz, V_inf, R, ...
    flight_phase, M_dot_FOD, M_dot_comp, solve_T, tol_T);

n_groups = numel(use_DNS_list) * numel(N_cool_list);
col      = lines(max(n_groups, 1));
g        = 0;
results_D = [];   % last discretised result, used for d_h_air below

for dns = use_DNS_list(:)'
    for N_cool_seg = N_cool_list(:)'
        g = g + 1;
        for i = 1:numel(N_list)
            N_segments = N_list(i);
            ls = line_styles{mod(i-1, numel(line_styles)) + 1};

            % The id is looked up fresh each time, because the first
            % variant that gets saved is what creates it.
            id_D = find_matching_id(inputs_D);

            res = [];
            from_cache = false;
            if ~isempty(id_D) && ~force_recompute
                candidate_file = fullfile('results', sprintf('id-%d', id_D), ...
                    disc_results_filename(N_segments, N_cool_seg, dns, id_D, solve_T));
                if exist(candidate_file, 'file')
                    fprintf('Found cached results, loading: %s\n', candidate_file);
                    loaded     = load(candidate_file, 'results');
                    res        = loaded.results;
                    from_cache = true;
                end
            end

            if isempty(res)
                res = run_disc_model_fwdpass(dns, N_segments, N_cool_seg, e, r, hx_theta, fan, fpr_init, ...
                    M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
                    d2_init, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, ...
                    solve_T, save_results, overwrite, tol_T);
                if dns, ran_dns_fresh = true; end
            end
            results_D = res;

            % Plot type / label
            if N_cool_seg == 0
                model_type = ternary(dns, 'dns', 'discretised');
                lbl = sprintf('%s 1D, N=%d (L=%.3fm)', ternary(dns, 'DNS', 'Corr.'), N_segments, res.L_solution);
            else
                model_type = ternary(dns, '2DNS', '2D');
                lbl = sprintf('%s 2D, N=%dx%d (L=%.3fm)', ternary(dns, 'DNS', 'Corr.'), N_segments, N_cool_seg, res.L_solution);
            end

            figs = plot_stations(res, model_type, 'figs', figs, 'color', col(g,:), ...
                'linestyle', ls, 'label', lbl);
            figs = plot_HX(res, model_type, 'figs', figs, 'color', col(g,:), ...
                'linestyle', ls, 'label', lbl);

            summary(end+1,:) = {lbl, res.L_solution, res.dp_hx, res.T4, res.M_hx, res.drag_HX, from_cache}; %#ok<SAGROW>
        end
    end
end

%%%--- SUMMARY TABLE ---%%%
fprintf("\n========================================\n");
fprintf("           MODEL COMPARISON             \n");
fprintf("========================================\n");
fprintf("%-34s %9s %10s %9s %9s %9s %7s\n", "Model", "L [m]", "dp_hx [Pa]", "T4 [K]", "M_hx [kg]", "drag [N]", "cached");
for k = 1:size(summary, 1)
    fprintf("%-34s %9.4f %10.2f %9.2f %9.2f %9.2f %7s\n", summary{k,1}, summary{k,2}, summary{k,3}, ...
        summary{k,4}, summary{k,5}, summary{k,6}, ternary(summary{k,7}, 'yes', 'no'));
end

%%%--- DNS COVERAGE (only meaningful for fresh DNS runs) ---%%%
if any_dns
    if ran_dns_fresh
        plot_thermoturb_coverage();
    else
        fprintf('Skipping DNS coverage plot: all DNS results were loaded from cache, so no DNS queries were logged.\n');
    end
end

%%%--- CFD (optional) ---%%%
if include_cfd && ~isempty(results_D)
    results_cfd = process_wall_CFD_data(...
        "C:\Users\tiesh\TUDelft\Thesis\CFD\\4by4_150ch_wall_q_air", ...
        "C:\Users\tiesh\TUDelft\Thesis\CFD\\4by4_150ch_wall_Twall_air", ...
        "C:\Users\tiesh\TUDelft\Thesis\CFD\\4by4_150ch_wall_Tbulk_air", ...
        results_D.d_h_air);
    % figs = plot_HX(results_cfd, 'cfd', 'figs', figs);
end

%% ---- Local helper: inline if ---- %%
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

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