%% compare_models.m  (model only)
% Compares the HX models themselves, on given HX inlet conditions, WITHOUT the
% rest of the TMS solve (no propeller, diffuser, fan, nozzle or pressure-balance
% iteration). The solved version of this script is compare_models_solved.m.
%
% What is called: HX_deltaT, then HX_design1 (lumped) and/or HX_design1_disc
% (discretised). Nothing else from the TMS chain.
%
% Inputs are the conditions at the HX inlet (station 3): mass flow, velocity,
% static pressure and temperature, plus the usual e and r and the coolant side.
% The HX frontal diameter follows from continuity: A3 = M_dot_3/(rho3*v3).
% To take these values from an earlier solved run, set solved_file.
%
% Choose what to compare, as in the solved script:
%   include_lumped  - add the lumped model
%   use_DNS_list    - correlation-based (false) and/or DNS-based (true)
%   N_cool_list     - 0 = old iterative 1D, >= 1 = 2D direct marching
%   N_list          - air-side segments
% Every combination of use_DNS_list x N_cool_list x N_list is one line in the figures.
%
% Saving / caching works as in the solved script (look for a saved result first,
% run only if none exists, save per save_results / overwrite, force_recompute),
% with these differences:
%   - file names start with 'model_only_', e.g. model_only_a50_3.mat
%   - they live in their own folder: results/model_only/id-<id>/
%   - the ids come from their own index, results/model_only/model_only_results_index.mat,
%     so the index used by the solved runs is not touched
%   - optionally a .json copy is written next to each .mat (save_json)
% plot_stations is not used: it needs the stations around the HX.

clear all
clc

%%%--- INPUT CASE: HX INLET (STATION 3) ---%%%

M_dot_3 = 5.0;       % air mass flow into the HX [kg/s]            (REPLACE with your case)
v3      = 60;        % air velocity at the HX inlet [m/s]          (REPLACE)
P3      = 60000;     % static pressure at the HX inlet [Pa]        (REPLACE)
T3      = 290;       % static temperature at the HX inlet [K]      (REPLACE)
P3_0    = [];        % total pressure [Pa]; [] = from P3, T3, v3 (isentropic)

% Optional: take M_dot_3, v3, P3, T3 (and P3_0) from a saved SOLVED result, e.g.
%   solved_file = fullfile('results', 'id-2', 'a50_2.mat');
solved_file = 'results/id-4/a50-c20_4.mat';

% Coolant side and heat duty
M_dot_coolant = 44.4;      % Max: 44.4, Min: 29
T_in_fc       = 70+273;
T_out_fc      = 85+273;
Q_tot         = 2.25e6;
n_modules     = 2;         % ducts per nacelle

% HX geometry
e        = 4 - 2*0.1*2/sqrt(3);
r        = 4 - 2*0.1*2/sqrt(3);
hx_theta = 60;

% Other
h        = 4876;     % altitude [m], passed to the lumped HX file
R        = 287;
solve_T  = 0;
tol_T    = 0;

%%%--- WHAT TO COMPARE ---%%%

include_lumped = true;
use_DNS_list   = [false true];
N_cool_list    = [0 5 20];
N_list         = 50;

line_styles = {'-', '--', ':', '-.'};   % one per entry of N_list, in order

%%%--- RUN / CACHE OPTIONS ---%%%

save_results    = true;
overwrite       = false;   % keep existing cached files; save new versions instead of clobbering
force_recompute = false;   % true = ignore saved results and re-run everything
save_json       = true;    % also write a text (.json) copy next to each saved .mat (opens in Notepad++)

include_cfd = true;       % true = also compare with the CFD meshes below

cfd_root = 'C:\Users\tiesh\TUDelft\Thesis\CFD\4by4_150ch_wall\run_2024_1';
% One row per mesh, COARSEST first: {legend label, mesh folder, run folder}
cfd_cases = {
'baseline', 'baseline',  'solution2/reports';
'fine',     'fine',      'solution2/reports';
'finest',   'fine_v2',   'solution2/reports';
};
cfd_export_files = {'q_air.srp', 'Twall_air.srp', 'Tbulk_air.srp', 'Pbulk_air.srp', 'Vbulk_air.srp', 'm_dot_air.srp'};
cfd_opt_names    = {'rho', 'mu', 'tau', 'yplus'};
cfd_opt_files    = {'rho_air.srp', 'mu_air.srp', 'tau_wall_air.srp', 'yplus_air.srp'};
include_coolant_cfd = false;   % also the coolant side (needs *_coolant.srp and Tbulk_coolant_in/_out.srp)
cfd_P_operating     = 55000;   % CFD operating pressure [Pa] (Fluent Pbulk is gauge); empty = shift to model inlet pressure
cfd_l_index         = 1;       % coolant-side column compared for the 2D models

%%%--- DERIVED HX INLET STATE ---%%%

if ~isempty(solved_file)
    S = load(solved_file, 'results');
    M_dot_3 = S.results.M_dot_3;  v3 = S.results.v3;
    P3 = S.results.P3;            T3 = S.results.T3;
    if isfield(S.results, 'P3_0'), P3_0 = S.results.P3_0; end
    fprintf('HX inlet state taken from %s\n', solved_file);
end

rho3 = rho_air(P3, T3);
A3   = M_dot_3/(rho3*v3);                 % continuity
d3   = sqrt(4*A3/pi);                     % diameter of the circle the HX frontal area is based on
cp3  = cp_air(T3);
T3_0 = T3 + v3^2/(2*cp3);
if isempty(P3_0)
    gam  = gamma_air(cp3);
    P3_0 = P3*(T3_0/T3)^(gam/(gam-1));
end
fprintf('HX inlet: M_dot_3 = %.4f kg/s, v3 = %.2f m/s, P3 = %.0f Pa, T3 = %.2f K, rho3 = %.4f kg/m3\n', ...
    M_dot_3, v3, P3, T3, rho3);
fprintf('          A3 = %.5f m2, d3 = %.4f m, P3_0 = %.0f Pa\n\n', A3, d3, P3_0);

% Everything the HX files need, in one struct (also the saving key)
in = struct('e', e, 'r', r, 'hx_theta', hx_theta, 'n_modules', n_modules, ...
    'M_dot_coolant', M_dot_coolant, 'Q_tot', Q_tot, 'T_in_fc', T_in_fc, 'T_out_fc', T_out_fc, ...
    'h', h, 'R', R, 'M_dot_3', M_dot_3, 'v3', v3, 'P3', P3, 'T3', T3, 'P3_0', P3_0, ...
    'solve_T', solve_T, 'tol_T', tol_T);
s3 = struct('A3', A3, 'd3', d3, 'T3_0', T3_0);

figs = struct();
summary = {};   % rows: label, L, dp_hx, T4, M_hx, drag_HX, from_cache

any_dns = any(use_DNS_list);
if any_dns
    thermoturb_query_log('reset');   % so the coverage plot only shows this script's fresh runs
end
ran_dns_fresh = false;

%%%--- LUMPED MODEL ---%%%
% Lumped has no tol_T: the lookup ignores it, so lumped shares the id of the
% discretised runs of the same inputs.

if include_lumped
    id_L = find_model_only_id(in, true);

    results_L = [];
    from_cache = false;
    if ~isempty(id_L) && ~force_recompute
        candidate_file = fullfile('results', 'model_only', sprintf('id-%d', id_L), ...
            ['model_only_' lumped_results_filename(id_L, solve_T)]);
        if exist(candidate_file, 'file')
            fprintf('Found cached results, loading: %s\n', candidate_file);
            loaded     = load(candidate_file, 'results');
            results_L  = loaded.results;
            from_cache = true;
        end
    end

    if isempty(results_L)
        results_L = run_model_only_lumped(in, s3);
        if save_results
            save_model_only(results_L, in, true, @(id) lumped_results_filename(id, solve_T), overwrite, save_json);
        end
    end

    figs = plot_HX(results_L, 'lumped', 'figs', figs);
    summary(end+1,:) = {'Lumped', results_L.L_solution, results_L.dp_hx, results_L.T4, ...
                        results_L.M_hx, results_L.drag_HX, from_cache};
end

%%%--- DISCRETISED VARIANTS ---%%%

n_groups = numel(use_DNS_list) * numel(N_cool_list);
col      = lines(max(n_groups, 1));
g        = 0;
results_D = [];   % last discretised result, used for d_h_air below
models_D  = {};   % every discretised result plotted, with its label (for the CFD comparison)
labels_D  = {};

for dns = use_DNS_list(:)'
    for N_cool_seg = N_cool_list(:)'
        g = g + 1;
        for i = 1:numel(N_list)
            N_segments = N_list(i);
            ls = line_styles{mod(i-1, numel(line_styles)) + 1};

            % The id is looked up fresh each time, because the first
            % variant that gets saved is what creates it.
            id_D = find_model_only_id(in, false);

            res = [];
            from_cache = false;
            if ~isempty(id_D) && ~force_recompute
                candidate_file = fullfile('results', 'model_only', sprintf('id-%d', id_D), ...
                    ['model_only_' disc_results_filename(N_segments, N_cool_seg, dns, id_D, solve_T)]);
                if exist(candidate_file, 'file')
                    fprintf('Found cached results, loading: %s\n', candidate_file);
                    loaded     = load(candidate_file, 'results');
                    res        = loaded.results;
                    from_cache = true;
                end
            end

            if isempty(res)
                res = run_model_only_disc(dns, N_segments, N_cool_seg, in, s3);
                if dns, ran_dns_fresh = true; end
                if save_results
                    save_model_only(res, in, false, ...
                        @(id) disc_results_filename(N_segments, N_cool_seg, dns, id, solve_T), overwrite, save_json);
                end
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

            figs = plot_HX(res, model_type, 'figs', figs, 'color', col(g,:), ...
                'linestyle', ls, 'label', lbl);

            summary(end+1,:) = {lbl, res.L_solution, res.dp_hx, res.T4, res.M_hx, res.drag_HX, from_cache}; %#ok<SAGROW>
            models_D{end+1} = res;  labels_D{end+1} = lbl; %#ok<SAGROW>
        end
    end
end

%%%--- SUMMARY TABLE ---%%%
fprintf("\n========================================\n");
fprintf("        MODEL-ONLY COMPARISON           \n");
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
if include_cfd && ~isempty(models_D)
    d_h = results_D.d_h_air;
    [cfd, cfd_cool, cfd_labels, cool_ends] = load_cfd_exports(cfd_root, cfd_cases, cfd_export_files, ...
        cfd_opt_names, cfd_opt_files, d_h, include_coolant_cfd);
    args = {'cfd_labels', cfd_labels, 'results_2D', models_D, 'labels_2D', labels_D, ...
            'cfd_coolant', cfd_cool, 'cfd_coolant_ends', cool_ends, 'cfd_P_operating', cfd_P_operating, 'l_index', cfd_l_index};
    if include_lumped, args = [args, {'lumped_results', results_L}]; end
    comparison = compare_cfd_meshes_vs_2d(cfd, args{:}); %#ok<NASGU>
end

%% ---- Local helper: inline if ---- %%
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

%% ---- Local helper: lumped HX model on the given HX inlet state ---- %%
function results = run_model_only_lumped(in, s3)
t_start = tic;
counter = 1;
x       = 1;

[T_mean_h, T_mean_c, T_c_i, T_c_o, T_h_i, T_h_o, C_h, C_c, C_star] = ...
    HX_deltaT(in.T3, in.M_dot_3, in.M_dot_coolant, in.T_out_fc, in.T_in_fc);

Pr_coolant = mu_EG(T_mean_h)*cp_EG_50_50(T_mean_h)/k_EG(T_mean_h);
Pr_air     = mu_air(T_mean_c)*cp_air(T_mean_c)/k_air(T_mean_c);

[dp_coolant, d_h_air, M_dot_4, ~, ~, dp_hx, N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, ...
    ~, ~, v_channel_air, ~, ~, ~, ~, Re_air, ~, h_air, ~, L_solution, ~, v4, P4_0, M4, T4, T4_0, P4, ...
    ~, M_hx, ~, drag_HX, Q_pred_sol, T_cool_out, Pr_air, f_air, Nu_air, f_hx, inlet_dp] = ...
    HX_design1(in.e, in.r, in.hx_theta, counter, x, s3.A3, in.v3, in.R, in.P3, in.P3_0, in.h, s3.d3, ...
    T_h_o, in.n_modules, T_h_i, T_c_i, T_c_o, T_mean_h, T_mean_c, in.Q_tot, C_h, C_c, C_star, ...
    Pr_air, Pr_coolant, in.M_dot_coolant, in.M_dot_3, in.T3);

results.L_solution      = L_solution;
results.dp_hx           = dp_hx;
results.dp_coolant      = dp_coolant;
results.T4              = T4;
results.M_hx            = M_hx;
results.drag_HX         = drag_HX;
results.T_cool_out      = T_cool_out;
results.P3              = in.P3;   results.P3_0 = in.P3_0;   results.P4 = P4;   results.P4_0 = P4_0;
results.T3              = in.T3;   results.T3_0 = s3.T3_0;   results.T4_0 = T4_0;
results.v3              = in.v3;   results.v4 = v4;
results.M4              = M4;
results.M_dot_3         = in.M_dot_3;
results.M_dot_4         = M_dot_4;
results.A3              = s3.A3;   results.d3 = s3.d3;
results.T_c_i           = T_c_i;
results.T_h_i           = T_h_i;
results.T_h_o           = T_h_o;
results.T_mean_c        = T_mean_c;
results.T_mean_h        = T_mean_h;
results.N_segments      = 1;
results.solve_T         = in.solve_T;
results.model_only      = true;
results.Q_pred_solution = Q_pred_sol;
results.h_air           = h_air;
results.v_channel_air   = v_channel_air;
results.f_air           = f_air;
results.Re_air          = Re_air;
results.Pr_air          = Pr_air;
results.Nu_air          = Nu_air;
results.f_hx            = f_hx;
results.inlet_dp        = inlet_dp;
results.d_h_air         = d_h_air;
results.N_fin_air       = N_fin_air;       results.N_fin_coolant = N_fin_coolant;
results.N_air_pass      = N_air_pass;      results.N_coolant_pass = N_coolant_pass;
results.t_compute       = toc(t_start);
fprintf('Lumped model: L = %.4f m  (%.3f s)\n', L_solution, results.t_compute);
end

%% ---- Local helper: discretised HX model on the given HX inlet state ---- %%
function results = run_model_only_disc(use_DNS, N_segments, N_cool_seg, in, s3)
t_start = tic;
counter = 1;

if use_DNS
    thermoturb_query_log('new_iter');
end

[T_mean_h, ~, T_c_i, ~, T_h_i, T_h_o, ~, ~, ~] = ...
    HX_deltaT(in.T3, in.M_dot_3, in.M_dot_coolant, in.T_out_fc, in.T_in_fc);

[dp_coolant_loop, d_h_air, M_dot_4, ~, ~, dp_hx, ...
    N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, ~, ~, ...
    ~, v_channel_coolant, d_h_coolant, ~, A_o_air, ...
    ~, Re_coolant, ~, ~, L_solution, ~, v4, P4_0, ...
    M4, T4, T4_0, P4, ~, M_hx, ~, drag_HX, T_cool_seg, dp_cool_seg, ...
    T_air_seg, P_air_seg, v_air_seg, Re_air_seg, Pr_air_seg, v_channel_seg, ...
    ~, f_air_seg, Nu_air_seg, h_air_seg, ~, ...
    ~, ~, ~, Q_seg_arr, dp_seg_arr, ...
    Q_pred_solution, T_h_o_solution, P_0_air_seg, T_0_air_seg, M_air_seg, f_hx_seg, inlet_dp, outlet_dp, ...
    delta_BL_seg, A_free_seg, d_h_bulk_seg, T_mean_c_arr, T_mean_h_arr, ...
    h_cool_seg] = ...
    HX_design1_disc(use_DNS, in.e, in.r, in.hx_theta, counter, s3.A3, in.v3, in.R, in.P3, ...
    s3.d3, T_h_o, in.n_modules, T_h_i, T_c_i, T_mean_h, ...
    in.Q_tot, in.M_dot_coolant, in.M_dot_3, in.T3, N_segments, N_cool_seg, in.tol_T);

% Air mass flow at every node, as in run_disc_model_fwdpass
m_dot_seg = zeros(1, N_segments+1);
for k = 1:N_segments+1
    rho_k        = rho_air(P_air_seg(k), T_air_seg(k));
    m_dot_seg(k) = rho_k * v_channel_seg(k) * A_o_air * N_fin_air * N_air_pass;
end
m_dot_seg(1) = m_dot_seg(2);

results.L_solution      = L_solution;
results.dp_hx           = dp_hx;
results.dp_coolant      = dp_coolant_loop;
results.T4              = T4;
results.M_hx            = M_hx;
results.drag_HX         = drag_HX;
results.T_cool_out      = T_h_o_solution;
results.P3              = in.P3;   results.P3_0 = in.P3_0;   results.P4 = P4;   results.P4_0 = P4_0;
results.T3              = in.T3;   results.T3_0 = s3.T3_0;   results.T4_0 = T4_0;
results.v3              = in.v3;   results.v4 = v4;
results.M4              = M4;
results.M_dot_3         = in.M_dot_3;
results.M_dot_4         = M_dot_4;
results.m_dot_seg       = m_dot_seg;
results.A3              = s3.A3;   results.d3 = s3.d3;
results.T_c_i           = T_c_i;
results.T_h_i           = T_h_i;
results.T_h_o           = T_h_o;
results.N_segments      = N_segments;
results.N_cool_seg      = N_cool_seg;
results.b_hx            = sqrt(pi*s3.d3^2/(4*cosd(in.hx_theta)));
results.use_DNS         = use_DNS;
results.solve_T         = in.solve_T;
results.model_only      = true;
results.T_air_seg       = T_air_seg;
results.P_air_seg       = P_air_seg;
results.v_air_seg       = v_air_seg;
results.v_channel_seg   = v_channel_seg;
results.Re_air_seg      = Re_air_seg;
results.Pr_air_seg      = Pr_air_seg;
results.f_air_seg       = f_air_seg;
results.Nu_air_seg      = Nu_air_seg;
results.h_air_seg       = h_air_seg;
results.h_cool_seg      = h_cool_seg;
results.v_cool_seg      = v_channel_coolant;
results.Re_cool_seg     = Re_coolant;
results.d_h_coolant     = d_h_coolant;
results.Q_seg_arr       = Q_seg_arr;
results.dp_seg_arr      = dp_seg_arr;
results.T_cool_seg      = T_cool_seg;
results.dp_cool_seg     = dp_cool_seg;
results.Q_pred_solution = Q_pred_solution;
results.P_0_air_seg     = P_0_air_seg;
results.T_0_air_seg     = T_0_air_seg;
results.M_air_seg       = M_air_seg;
results.f_hx_seg        = f_hx_seg;
results.inlet_dp        = inlet_dp;
results.outlet_dp       = outlet_dp;
results.d_h_air         = d_h_air;
results.delta_BL_seg    = delta_BL_seg;
results.A_free_seg      = A_free_seg;
results.d_h_bulk_seg    = d_h_bulk_seg;
results.T_mean_h_arr    = T_mean_h_arr;
results.T_mean_c_arr    = T_mean_c_arr;
results.A_o_air         = A_o_air;
results.N_fin_air       = N_fin_air;       results.N_fin_coolant = N_fin_coolant;
results.N_air_pass      = N_air_pass;      results.N_coolant_pass = N_coolant_pass;
results.t_compute       = toc(t_start);
if use_DNS
    results.thermoturb_log = thermoturb_query_log('get');
end
fprintf('Discretised model (DNS=%d, N=%d, N_cool=%d): L = %.4f m  (%.3f s)\n', ...
    use_DNS, N_segments, N_cool_seg, L_solution, results.t_compute);
end

%% ---- Local helper: look up the id of an input case in the model-only index ---- %%
% ignore_tol_T = true (lumped): rows are compared with tol_T removed on both sides.
function id = find_model_only_id(in, ignore_tol_T)
index_file = fullfile('results', 'model_only', 'model_only_results_index.mat');
id = [];
if ~exist(index_file, 'file'), return, end
loaded = load(index_file, 'model_only_index');
idx = loaded.model_only_index;
for row = 1:height(idx)
    si = idx.inputs{row};
    if ignore_tol_T
        match = isequal(strip_tol(si), strip_tol(in));
    else
        match = isequal(si, in);
    end
    if match
        id = idx.id(row);
        return
    end
end
end

%% ---- Local helper: save one model-only result (same scheme as the solved runs) ---- %%
% name_fn(id) returns the plain results filename (disc_results_filename or
% lumped_results_filename); 'model_only_' goes in front of it and the file is put
% in results/model_only/id-<id>/. With overwrite = false an existing file is kept
% and the new one gets a -1, -2, ... version suffix (next_versioned_filename).
function results_file = save_model_only(results, in, is_lumped, name_fn, overwrite, save_json)
results_dir = fullfile('results', 'model_only');
index_file  = fullfile(results_dir, 'model_only_results_index.mat');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

if exist(index_file, 'file')
    loaded = load(index_file, 'model_only_index');
    idx = loaded.model_only_index;
else
    idx = table('Size', [0 2], 'VariableTypes', {'double', 'cell'}, 'VariableNames', {'id', 'inputs'});
end

% Same inputs -> same id. The lumped model has no tol_T, so it matches a
% discretised row of the same inputs; a lumped-only row is upgraded to the
% full key when a discretised run is saved for it.
id = [];
for row = 1:height(idx)
    si = idx.inputs{row};
    if is_lumped
        match = isequal(strip_tol(si), strip_tol(in));
    else
        match = isequal(si, in);
    end
    if match
        id = idx.id(row);
        break
    end
end
if isempty(id) && ~is_lumped
    for row = 1:height(idx)
        si = idx.inputs{row};
        if ~isfield(si, 'tol_T') && isequal(si, strip_tol(in))
            id = idx.id(row);
            idx.inputs{row} = in;
            model_only_index = idx;
            save(index_file, 'model_only_index');
            break
        end
    end
end
if isempty(id)
    id = height(idx) + 1;
    idx = [idx; table(id, {in}, 'VariableNames', {'id', 'inputs'})];
    model_only_index = idx;
    save(index_file, 'model_only_index');
end

id_dir = fullfile(results_dir, sprintf('id-%d', id));
if ~exist(id_dir, 'dir'), mkdir(id_dir); end

base_fname = ['model_only_' name_fn(id)];
if overwrite
    results_file = fullfile(id_dir, base_fname);
else
    results_file = fullfile(id_dir, next_versioned_filename(id_dir, base_fname));
end

save(results_file, 'results');
fprintf('Saved results: %s\n', results_file);

if save_json
    try
        r2 = results;
        if isfield(r2, 'thermoturb_log'), r2 = rmfield(r2, 'thermoturb_log'); end
        try
            txt = jsonencode(r2, 'PrettyPrint', true);   % R2021a and later
        catch
            txt = jsonencode(r2);
        end
        json_file = [results_file(1:end-4) '.json'];
        fid = fopen(json_file, 'w');
        fwrite(fid, txt, 'char');
        fclose(fid);
        fprintf('Saved text copy: %s\n', json_file);
    catch ME
        warning('compare_models:json', 'Could not write the .json copy: %s', ME.message);
    end
end
end

%% ---- Local helper: remove tol_T from an input struct ---- %%
function s = strip_tol(s)
if isfield(s, 'tol_T'), s = rmfield(s, 'tol_T'); end
end