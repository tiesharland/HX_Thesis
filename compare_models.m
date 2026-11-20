%% compare_models.m  (model only)
% Compares the HX models themselves, on given HX inlet conditions, WITHOUT the
% rest of the TMS solve (no propeller, diffuser, fan, nozzle or pressure-balance
% iteration). The solved version of this script is compare_models_solved.m.
%
% What is called: HX_deltaT, then one or more lumped HX files (HX_design1, and the
% variants listed in lumped_variants, e.g. HX_design1_v2) and/or
% HX_design1_disc (discretised). Nothing else from the TMS chain.
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
% Load or recompute: set 'recompute' in the run options.
%   'none'            - load every model that has a saved result (default)
%   'all'             - recompute every model
%   {'v2', 'DNS&2D'}  - recompute only the models whose legend name contains one of these
%                       (not case sensitive; 'a&b' = both words must be in the name).
%                       All other models are loaded.
% A model with no saved result is always computed.
%
% Saving / caching works as in the solved script (look for a saved result first,
% run only if none exists, save per save_results / overwrite),
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
solved_file = 'results/id-4/a50_4.mat';

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
s_bend_switch = "NO";  % only passed to lumped files that take it (e.g. the older HX_design1 versions with an s-bend switch)
R        = 287;
solve_T  = 0;
tol_T    = 0;

%%%--- WHAT TO COMPARE ---%%%

include_lumped = true;

% Lumped HX files to compare. One row each, any function with the HX_design1 interface
% (28 inputs, or 29 with a leading s_bend_switch; 44 outputs (older files: 41, 40 or 33) -- see call_lumped_hx).
%   {function name, file tag, legend name, line style, colour}
% The file tag goes into the saved file name (lumped-<tag>_<id>.mat); '' = plain lumped_<id>.mat.
% The legend name is also what 'recompute' (run options) is matched against.
lumped_variants = {
    'HX_design1',       '',     'Lumped',        '-',  [0 0 0];
    'HX_design1_v2',    'v2',   'Lumped v2',     '--', [0.55 0.27 0.07];
};
lumped_use = [1 2];   % rows of lumped_variants to run; a row whose file is not on the path is skipped with a warning
use_DNS_list   = [false];
N_cool_list    = [0 20 318];
N_list         = 50;

% Discretised HX files (same interface as HX_design1_disc, incl. use_DNS). One row each:
%   {function name, file tag}
% The tag is added to the legend ('test') and the saved file name (a50-c5_dns-test_<id>.mat).
disc_variants = {
    'HX_design1_disc',       '';
    'HX_design1_disc_test',  'test';
};
disc_use = [1 2];     % rows of disc_variants to run; a row whose file is not on the path is skipped with a warning

line_styles = {'-', '--', ':', '-.'};   % one per entry of N_list, in order

%%%--- RUN / CACHE OPTIONS ---%%%

save_results    = false;
overwrite       = false;   % keep existing cached files; save new versions instead of clobbering

% Which models to load from saved results and which to recompute:
%   'none'  = load all that exist          'all' = recompute all
%   {'test'} or {'v2','DNS&2D'} = recompute only models whose legend name contains
%   one of these words ('&' = all words); everything else is loaded.
% Names look like: 'Lumped', 'Lumped v2', 'Corr. 1D, N=50', 'DNS 2D test, N=50x5'.
% Note: with overwrite = false a recomputed result is saved as a new version (-1, -2, ...)
% and the next run still loads the old file; set overwrite = true to replace it.
recompute       = 'all';

% Which saved version to load (only used when a model is loaded, not recomputed).
% Versions are the files with a -1, -2, ... suffix, e.g. model_only_a50_3-1.mat.
%   load_version = 0         the original file (no suffix)
%   load_version = 2         the file ending -2
%   load_version = 'latest'  the highest version that exists
% load_version_for overrides this for some models; same words as 'recompute',
% first matching row wins:
%   load_version_for = {'test', 1; 'v2', 'latest'};
% If the chosen version does not exist, a warning is shown and the highest older version
% (or the original file) is loaded; the model is only computed if nothing is saved.
load_version     = 'latest';
load_version_for = {};
save_json       = true;    % also write a text (.json) copy next to each saved .mat (opens in Notepad++)
show_plots      = false;    % false = skip plot_HX (and the DNS coverage plot): just run/load, save and print the summary
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
include_coolant_cfd = true;   % also the coolant side (needs *_coolant.srp and Tbulk_coolant_in/_out.srp)
cfd_P_operating     = 55000;   % CFD operating pressure [Pa] (Fluent Pbulk is gauge); empty = shift to model inlet pressure
cfd_l_index         = 1;       % coolant-side column compared for the 2D models
cfd_nu_ratio        = false;   % true = also the Nu_CFD/Nu_developed figure; its DNS part queries the thermoturb table

%%%--- SKIP VARIANT FILES THAT ARE NOT SUPPLIED ---%%%
lumped_use = existing_variants(lumped_variants, lumped_use);
disc_use   = existing_variants(disc_variants, disc_use);

recompute = cellstr(recompute);            % 'all' / 'none' / {patterns}
rc_hits   = false(1, numel(recompute));    % which patterns matched a model (warn if none)
if isempty(recompute) || strcmpi(recompute{1}, 'none')
    fprintf('Recompute: none (saved results are loaded where they exist)\n');
elseif strcmpi(recompute{1}, 'all')
    fprintf('Recompute: all models\n');
else
    fprintf('Recompute: models matching %s\n', strjoin(recompute, ', '));
end

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
ran_dns_fresh = false;

%%%--- LUMPED MODEL ---%%%
% Lumped has no tol_T: the lookup ignores it, so lumped shares the id of the
% discretised runs of the same inputs.

if include_lumped
    results_L = {};     % every lumped variant that was run (for the CFD comparison)
    names_L   = {};
    for iv = lumped_use(:)'
        fn_L    = lumped_variants{iv,1};   tag_L   = lumped_variants{iv,2};
        name_L  = lumped_variants{iv,3};
        ls_L    = lumped_variants{iv,4};   col_L   = lumped_variants{iv,5};
        id_L = find_model_only_id(in, true);
        [do_rc, hit] = want_recompute(recompute, name_L);  rc_hits = rc_hits | hit;

        res_L = [];
        from_cache = false;
        if ~isempty(id_L) && ~do_rc
            candidate_file = pick_version(fullfile('results', 'model_only', sprintf('id-%d', id_L)), ...
                ['model_only_' lumped_name(id_L, solve_T, tag_L)], version_for(load_version_for, load_version, name_L));
            if ~isempty(candidate_file)
                fprintf('Found cached results, loading: %s\n', candidate_file);
                loaded     = load(candidate_file, 'results');
                res_L      = loaded.results;
                from_cache = true;
            end
        end

        if isempty(res_L)
            res_L = run_model_only_lumped(in, s3, fn_L, s_bend_switch);
            if save_results
                save_model_only(res_L, in, true, @(id) lumped_name(id, solve_T, tag_L), overwrite, save_json);
            end
        end
        results_L{end+1} = res_L;  names_L{end+1} = name_L; %#ok<SAGROW>

        lbl_L = sprintf('%s (L=%.3fm)', name_L, res_L.L_solution);
        if show_plots
            figs = plot_HX(res_L, 'lumped', 'figs', figs, 'color', col_L, 'linestyle', ls_L, 'label', lbl_L);
        end
        summary(end+1,:) = {name_L, res_L.L_solution, res_L.dp_hx, res_L.T4, ...
                            res_L.M_hx, res_L.drag_HX, from_cache}; %#ok<SAGROW>
    end
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
          for iv = disc_use(:)'
            N_segments = N_list(i);
            ls = line_styles{mod(i-1, numel(line_styles)) + 1};
            fn_D = disc_variants{iv,1};  tag_D = disc_variants{iv,2};

            % Plot type / label (the variant tag, e.g. 'test', is added to the label only)
            tag_txt = ternary(isempty(tag_D), '', [' ' tag_D]);
            if N_cool_seg == 0
                model_type = ternary(dns, 'dns', 'discretised');
                lbl_base = sprintf('%s 1D%s, N=%d', ternary(dns, 'DNS', 'Corr.'), tag_txt, N_segments);
            else
                model_type = ternary(dns, '2DNS', '2D');
                lbl_base = sprintf('%s 2D%s, N=%dx%d', ternary(dns, 'DNS', 'Corr.'), tag_txt, N_segments, N_cool_seg);
            end
            [do_rc, hit] = want_recompute(recompute, lbl_base);  rc_hits = rc_hits | hit;

            % The id is looked up fresh each time, because the first
            % variant that gets saved is what creates it.
            id_D = find_model_only_id(in, false);

            res = [];
            from_cache = false;
            if ~isempty(id_D) && ~do_rc
                candidate_file = pick_version(fullfile('results', 'model_only', sprintf('id-%d', id_D)), ...
                    ['model_only_' disc_name(N_segments, N_cool_seg, dns, id_D, solve_T, tag_D)], ...
                    version_for(load_version_for, load_version, lbl_base));
                if ~isempty(candidate_file)
                    fprintf('Found cached results, loading: %s\n', candidate_file);
                    loaded     = load(candidate_file, 'results');
                    res        = loaded.results;
                    from_cache = true;
                end
            end

            if isempty(res)
                if dns && ~ran_dns_fresh
                    thermoturb_query_log('reset');   % first fresh DNS run: coverage plot shows only this script's runs
                end
                res = run_model_only_disc(dns, N_segments, N_cool_seg, in, s3, fn_D);
                if dns, ran_dns_fresh = true; end
                if save_results
                    save_model_only(res, in, false, ...
                        @(id) disc_name(N_segments, N_cool_seg, dns, id, solve_T, tag_D), overwrite, save_json);
                end
            end
            results_D = res;

            lbl = sprintf('%s (L=%.3fm)', lbl_base, res.L_solution);

            if show_plots
                figs = plot_HX(res, model_type, 'figs', figs, 'color', col(g,:), ...
                    'linestyle', ls, 'label', lbl);
            end

            summary(end+1,:) = {lbl, res.L_solution, res.dp_hx, res.T4, res.M_hx, res.drag_HX, from_cache}; %#ok<SAGROW>
            models_D{end+1} = res;  labels_D{end+1} = lbl; %#ok<SAGROW>
          end
        end
    end
end

if ~isempty(recompute) && ~any(strcmpi(recompute{1}, {'all', 'none'}))
    for k = find(~rc_hits)
        warning('compare_models:recomputeNoMatch', 'recompute entry "%s" matched no model that was run.', recompute{k});
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
if any_dns && show_plots
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
            'cfd_coolant', cfd_cool, 'cfd_coolant_ends', cool_ends, 'cfd_P_operating', cfd_P_operating, 'l_index', cfd_l_index, 'nu_ratio', cfd_nu_ratio};
    if include_lumped && ~isempty(results_L), args = [args, {'lumped_results', results_L, 'lumped_labels', names_L}]; end
    comparison = compare_cfd_meshes_vs_2d(cfd, args{:}); %#ok<NASGU>
end

%% ---- Local helper: keep only the variant rows whose function file is on the path ---- %%
% Row 1 (the original HX_design1 / HX_design1_disc) is never skipped silently: if it is
% missing the later call fails as it normally would.
function use = existing_variants(variants, use)
keep = true(size(use));
for k = 1:numel(use)
    if use(k) == 1, continue, end
    if exist(variants{use(k),1}, 'file') == 0
        warning('compare_models:variantMissing', '%s not found on the path -- skipped.', variants{use(k),1});
        keep(k) = false;
    end
end
use = use(keep);
end

%% ---- Local helper: should this model be recomputed? ---- %%
% recompute: {'all'}, {'none'} or a list of patterns; label: the model's legend name
% (without the length). A pattern 'a&b' needs both words in the label. hit marks the
% patterns that matched, so unused ones can be reported.
function [tf, hit] = want_recompute(recompute, label)
hit = false(1, numel(recompute));
if isempty(recompute) || strcmpi(recompute{1}, 'none'), tf = false; return, end
if strcmpi(recompute{1}, 'all'), tf = true; return, end
for k = 1:numel(recompute)
    if isempty(strtrim(recompute{k})), continue, end
    parts = strtrim(strsplit(recompute{k}, '&'));
    hit(k) = all(contains(label, parts, 'IgnoreCase', true));
end
tf = any(hit);
end

%% ---- Local helper: which version to load for this model ---- %%
% rules: {pattern, version; ...} (patterns as in 'recompute'); default: load_version.
function v = version_for(rules, default, label)
v = default;
for k = 1:size(rules, 1)
    if want_recompute(cellstr(rules{k,1}), label)
        v = rules{k,2};
        return
    end
end
end

%% ---- Local helper: path of the requested saved version, '' if there is none ---- %%
% v = 0: base file; v = n: base-n; v = 'latest': highest existing (base if no suffixed one).
function file = pick_version(id_dir, base_fname, v)
[~, name, ext] = fileparts(base_fname);
file = '';
if (ischar(v) || isstring(v)) && strcmpi(v, 'latest')
    % highest existing version number, even if lower ones (e.g. -1) are missing
    n = 0;
    d = dir(fullfile(id_dir, [name '-*' ext]));
    for k = 1:numel(d)
        tok = regexp(d(k).name, ['^' regexptranslate('escape', name) '-(\d+)' regexptranslate('escape', ext) '$'], 'tokens', 'once');
        if ~isempty(tok), n = max(n, str2double(tok{1})); end
    end
    v = n;
end
if v == 0
    cand = fullfile(id_dir, base_fname);
else
    cand = fullfile(id_dir, sprintf('%s-%d%s', name, v, ext));
end
if exist(cand, 'file')
    file = cand;
elseif v ~= 0
    % Requested version is missing: use the highest older one (or the original file) instead
    % of recomputing the model, and say so.
    n = v - 1;
    while n > 0 && ~exist(fullfile(id_dir, sprintf('%s-%d%s', name, n, ext)), 'file'), n = n - 1; end
    if n > 0
        alt = fullfile(id_dir, sprintf('%s-%d%s', name, n, ext));
    else
        alt = fullfile(id_dir, base_fname);
    end
    if exist(alt, 'file')
        warning('compare_models:versionMissing', '%s not found -- loading %s instead.', cand, alt);
        file = alt;
    else
        warning('compare_models:versionMissing', '%s not found and nothing older saved -- computing instead.', cand);
    end
end
end

%% ---- Local helper: inline if ---- %%
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

%% ---- Local helper: lumped HX model on the given HX inlet state ---- %%
% fname: name of the lumped HX function (HX_design1, HX_design1_v2, ...).
function results = run_model_only_lumped(in, s3, fname, s_bend_switch)
t_start = tic;
counter = 1;
x       = 1;

[b.T_mean_h, b.T_mean_c, b.T_c_i, b.T_c_o, b.T_h_i, b.T_h_o, b.C_h, b.C_c, b.C_star] = ...
    HX_deltaT(in.T3, in.M_dot_3, in.M_dot_coolant, in.T_out_fc, in.T_in_fc);

Pr_coolant = mu_EG(b.T_mean_h)*cp_EG_50_50(b.T_mean_h)/k_EG(b.T_mean_h);
Pr_air     = mu_air(b.T_mean_c)*cp_air(b.T_mean_c)/k_air(b.T_mean_c);

o = call_lumped_hx(fname, in, s3, b, counter, x, Pr_air, Pr_coolant, s_bend_switch);

results.L_solution      = o.L_solution;
results.dp_hx           = o.dp_hx;
results.dp_coolant      = o.dp_coolant_loop;
results.T4              = o.T4;
results.M_hx            = o.M_hx;
results.drag_HX         = o.drag_HX;
results.T_cool_out      = o.T_cool_out;
results.P3              = in.P3;   results.P3_0 = in.P3_0;   results.P4 = o.P4;   results.P4_0 = o.P4_0;
results.T3              = in.T3;   results.T3_0 = s3.T3_0;   results.T4_0 = o.T4_0;
results.v3              = in.v3;   results.v4 = o.v4;
results.M4              = o.M4;
results.M_dot_3         = in.M_dot_3;
results.M_dot_4         = o.M_dot_4;
results.A3              = s3.A3;   results.d3 = s3.d3;
results.T_c_i           = b.T_c_i;
results.T_h_i           = b.T_h_i;
results.T_h_o           = b.T_h_o;
results.T_mean_c        = b.T_mean_c;
results.T_mean_h        = b.T_mean_h;
results.N_segments      = 1;
results.solve_T         = in.solve_T;
results.model_only      = true;
results.hx_function     = fname;
results.filled_fields   = o.filled;     % outputs this HX file does not return, filled by call_lumped_hx
results.Q_pred_solution = o.Q_pred_solution;
results.h_air           = o.h_air;
results.v_channel_air   = o.v_channel_air;
results.f_air           = o.f_air;
results.Re_air          = o.Re_air;
results.Pr_air          = o.Pr_air;
results.Nu_air          = o.Nu_air;
results.f_hx            = o.f_hx;
results.inlet_dp        = o.inlet_dp;
results.outlet_dp       = o.outlet_dp;
results.Pr_coolant      = o.Pr_coolant;
results.f_coolant       = o.f_coolant;       % Darcy, as used for the coolant pressure drop
results.Nu_coolant      = o.Nu_coolant;
results.Re_coolant      = o.Re_coolant;
results.h_coolant       = o.h_coolant;
results.v_channel_coolant = o.v_channel_coolant;
results.d_h_coolant     = o.d_h_coolant;
results.d_h_air         = o.d_h_air;
results.N_fin_air       = o.N_fin_air;       results.N_fin_coolant = o.N_fin_coolant;
results.N_air_pass      = o.N_air_pass;      results.N_coolant_pass = o.N_coolant_pass;
results.t_compute       = toc(t_start);
fprintf('Lumped model %s: L = %.4f m  (%.3f s)\n', fname, o.L_solution, results.t_compute);
end

%% ---- Local helper: call a lumped HX function whatever its interface (28/29 inputs, 44/41/40/33 outputs) ---- %%
% Inputs: the standard HX_design1 list; if the file takes one more input, it is the
% leading s_bend_switch. Outputs: the standard 44 in the order below; a file that
% returns fewer gets the rest filled in (listed in o.filled):
%   Q_pred_solution = Q_tot/n_modules   (the length is solved so the prediction equals this)
%   T_cool_out      = T_h_i - Q_pred_solution/C_h
%   Pr_air          = the input value;  f_air, Nu_air, f_hx, inlet_dp, outlet_dp, f_coolant, Nu_coolant = NaN;  Pr_coolant = the input value
function o = call_lumped_hx(fname, in, s3, b, counter, x, Pr_air, Pr_coolant, s_bend_switch)
args = {in.e, in.r, in.hx_theta, counter, x, s3.A3, in.v3, in.R, in.P3, in.P3_0, in.h, s3.d3, ...
    b.T_h_o, in.n_modules, b.T_h_i, b.T_c_i, b.T_c_o, b.T_mean_h, b.T_mean_c, in.Q_tot, ...
    b.C_h, b.C_c, b.C_star, Pr_air, Pr_coolant, in.M_dot_coolant, in.M_dot_3, in.T3};
n_in = nargin(fname);
if n_in == numel(args) + 1
    args = [{s_bend_switch}, args];
elseif n_in ~= numel(args)
    error('compare_models:lumpedInputs', '%s takes %d inputs; expected %d (or %d with a leading s_bend_switch).', ...
        fname, n_in, numel(args), numel(args)+1);
end
n_out = nargout(fname);
if n_out < 1
    error('compare_models:lumpedOutputs', 'Cannot determine the number of outputs of %s.', fname);
end
names = {'dp_coolant_loop', 'd_h_air', 'M_dot_4', 'b_t_air', 'b_t_coolant', 'dp_hx', 'N_fin_air', ...
    'N_fin_coolant', 'N_air_pass', 'N_coolant_pass', 'NTU', 'R_tot', 'v_channel_air', 'v_channel_coolant', ...
    'd_h_coolant', 'A_o_coolant', 'A_o_air', 'Re_air', 'Re_coolant', 'h_air', 'h_coolant', 'L_solution', ...
    'UA_unit', 'v4', 'P4_0', 'M4', 'T4', 'T4_0', 'P4', 'F_drag_hx', 'M_hx', 'A4', 'drag_HX', ...
    'Q_pred_solution', 'T_cool_out', 'Pr_air', 'f_air', 'Nu_air', 'f_hx'};
if n_out >= 44
    % newer files: coolant Pr, friction factor and Nu come before inlet_dp / outlet_dp
    names = [names, {'Pr_coolant', 'f_coolant', 'Nu_coolant', 'inlet_dp', 'outlet_dp'}];
else
    names = [names, {'inlet_dp', 'outlet_dp'}];
end
outs = cell(1, n_out);
[outs{:}] = feval(fname, args{:});
o = struct();
for i = 1:min(n_out, numel(names))
    o.(names{i}) = outs{i};
end
o.filled = {};
if ~isfield(o, 'Q_pred_solution'), o.Q_pred_solution = in.Q_tot/in.n_modules; o.filled{end+1} = 'Q_pred_solution'; end
if ~isfield(o, 'T_cool_out'),      o.T_cool_out = b.T_h_i - o.Q_pred_solution/b.C_h; o.filled{end+1} = 'T_cool_out'; end
if ~isfield(o, 'Pr_air'),          o.Pr_air = Pr_air;  o.filled{end+1} = 'Pr_air'; end
if ~isfield(o, 'Pr_coolant'),      o.Pr_coolant = Pr_coolant;  o.filled{end+1} = 'Pr_coolant'; end
for nm = {'f_air', 'Nu_air', 'f_hx', 'inlet_dp', 'outlet_dp', 'f_coolant', 'Nu_coolant'}
    if ~isfield(o, nm{1}), o.(nm{1}) = NaN; o.filled{end+1} = nm{1}; end %#ok<AGROW>
end
if ~isempty(o.filled)
    fprintf('%s does not return: %s -- filled in (see call_lumped_hx).\n', fname, strjoin(o.filled, ', '));
end
end

%% ---- Local helper: discretised results filename, with an optional variant tag ---- %%
% tag '' -> as disc_results_filename; tag 'test' -> a50-c5_dns-test_<id>.mat
function fname = disc_name(N_segments, N_cool_seg, use_DNS, id, solve_T, tag)
fname = disc_results_filename(N_segments, N_cool_seg, use_DNS, id, solve_T);
if ~isempty(tag)
    fname = regexprep(fname, '_(\d+)(-T)?\.mat$', ['-' tag '_$1$2.mat']);
end
end

%% ---- Local helper: lumped results filename, with an optional variant tag ---- %%
% tag '' -> lumped_<id>.mat (as lumped_results_filename); tag 'v2' -> lumped-v2_<id>.mat
function fname = lumped_name(id, solve_T, tag)
fname = lumped_results_filename(id, solve_T);
if ~isempty(tag)
    fname = strrep(fname, 'lumped_', ['lumped-' tag '_']);
end
end

%% ---- Local helper: discretised HX model on the given HX inlet state ---- %%
% fname: name of the discretised HX function (HX_design1_disc, HX_design1_disc_test, ...).
function results = run_model_only_disc(use_DNS, N_segments, N_cool_seg, in, s3, fname)
t_start = tic;
n_out = nargout(fname);
if nargin(fname) ~= 22 || n_out ~= 68
    error('compare_models:discInterface', ...
        '%s must have the interface of HX_design1_disc (22 inputs, 68 outputs; the outputs were renamed -- see HX_design1_disc.m); it has %d inputs and %d outputs.', ...
        fname, nargin(fname), n_out);
end
counter = 1;

if use_DNS
    thermoturb_query_log('new_iter');
end

[T_mean_h, ~, T_c_i, ~, T_h_i, T_h_o, ~, ~, ~] = ...
    HX_deltaT(in.T3, in.M_dot_3, in.M_dot_coolant, in.T_out_fc, in.T_in_fc);

% Outputs in the order of HX_design1_disc.
names = {'dp_coolant_loop', 'd_h_air', 'M_dot_4', 'b_t_air', 'b_t_coolant', 'dp_hx', ...
    'N_fin_air', 'N_fin_coolant', 'N_air_pass', 'N_coolant_pass', 'NTU', 'R_tot', ...
    'v_channel_air', 'v_cool_seg', 'd_h_coolant', 'A_o_coolant', 'A_o_air', ...
    'Re_air', 'Re_cool_seg', 'h_air', 'h_coolant', 'L_solution', 'UA_unit', 'v4', 'P4_0', ...
    'M4', 'T4', 'T4_0', 'P4', 'F_drag_hx', 'M_hx', 'A4', 'drag_HX', 'T_cool_seg', 'dp_cool_seg', ...
    'T_air_seg', 'P_air_seg', 'Re_air_seg', 'Pr_air_seg', 'v_air_seg', ...
    'K_seg', 'Cf_air_seg', 'Nu_air_seg', 'h_air_seg', 'eta_fin_seg', ...
    'UA_seg_arr', 'NTU_seg_arr', 'eps_seg_arr', 'Q_seg_arr', 'dp_seg_arr', ...
    'Q_pred_solution', 'T_h_o_solution', 'P_0_air_seg', 'T_0_air_seg', 'M_air_seg', 'f_air_seg', 'inlet_dp', 'outlet_dp', ...
    'delta_BL_seg', 'A_free_seg', 'd_h_bulk_seg', 'T_mean_c_arr', 'T_mean_h_arr', ...
    'h_cool_seg', 'Pr_cool_seg', 'Cf_coolant_seg', 'Nu_cool_seg', 'f_coolant_seg'};
outs = cell(1, n_out);
[outs{:}] = feval(fname, use_DNS, in.e, in.r, in.hx_theta, counter, s3.A3, in.v3, in.R, in.P3, ...
    s3.d3, T_h_o, in.n_modules, T_h_i, T_c_i, T_mean_h, ...
    in.Q_tot, in.M_dot_coolant, in.M_dot_3, in.T3, N_segments, N_cool_seg, in.tol_T);
o = cell2struct(outs(:), names(1:n_out)', 1);
dp_coolant_loop = o.dp_coolant_loop;  d_h_air = o.d_h_air;  M_dot_4 = o.M_dot_4;  dp_hx = o.dp_hx;
N_fin_air = o.N_fin_air;  N_fin_coolant = o.N_fin_coolant;  N_air_pass = o.N_air_pass;  N_coolant_pass = o.N_coolant_pass;
v_channel_coolant = o.v_cool_seg;  d_h_coolant = o.d_h_coolant;  A_o_air = o.A_o_air;
Re_coolant = o.Re_cool_seg;  L_solution = o.L_solution;  v4 = o.v4;  P4_0 = o.P4_0;  M4 = o.M4;
T4 = o.T4;  T4_0 = o.T4_0;  P4 = o.P4;  M_hx = o.M_hx;  drag_HX = o.drag_HX;
T_cool_seg = o.T_cool_seg;  dp_cool_seg = o.dp_cool_seg;  T_air_seg = o.T_air_seg;  P_air_seg = o.P_air_seg;
v_air_seg = o.v_air_seg;  Re_air_seg = o.Re_air_seg;  Pr_air_seg = o.Pr_air_seg;
Nu_air_seg = o.Nu_air_seg;  h_air_seg = o.h_air_seg;
Q_seg_arr = o.Q_seg_arr;  dp_seg_arr = o.dp_seg_arr;
Q_pred_solution = o.Q_pred_solution;  T_h_o_solution = o.T_h_o_solution;
P_0_air_seg = o.P_0_air_seg;  T_0_air_seg = o.T_0_air_seg;  M_air_seg = o.M_air_seg;
inlet_dp = o.inlet_dp;  outlet_dp = o.outlet_dp;
delta_BL_seg = o.delta_BL_seg;  A_free_seg = o.A_free_seg;  d_h_bulk_seg = o.d_h_bulk_seg;
T_mean_c_arr = o.T_mean_c_arr;  T_mean_h_arr = o.T_mean_h_arr;  h_cool_seg = o.h_cool_seg;

% Air mass flow at every node, as in run_disc_model_fwdpass
m_dot_seg = zeros(1, N_segments+1);
for k = 1:N_segments+1
    rho_k        = rho_air(P_air_seg(k), T_air_seg(k));
    m_dot_seg(k) = rho_k * v_air_seg(k) * A_o_air * N_fin_air * N_air_pass;
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
results.hx_function     = fname;
results.T_air_seg       = T_air_seg;
results.P_air_seg       = P_air_seg;
results.v_air_seg       = v_air_seg;
results.Re_air_seg      = Re_air_seg;
results.Pr_air_seg      = Pr_air_seg;
results.Cf_air_seg      = o.Cf_air_seg;      % Fanning
results.f_air_seg       = o.f_air_seg;       % Darcy
results.Nu_air_seg      = Nu_air_seg;
results.h_air_seg       = h_air_seg;
results.h_cool_seg      = h_cool_seg;
results.Pr_cool_seg     = o.Pr_cool_seg;
results.Cf_coolant_seg  = o.Cf_coolant_seg;  % Fanning
results.f_coolant_seg   = o.f_coolant_seg;   % Darcy, as used for dp_cool
results.Nu_cool_seg     = o.Nu_cool_seg;
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