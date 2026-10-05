%% Ncoolseg_L_curves.m
% How does the sized HX length L_solution change with N_cool_seg (the
% coolant-side resolution of the 2D direct-marching path)?
%
%   One COLOUR per input condition, one LINE STYLE per N_segments (solid =
%   finest). A second figure shows the same curves as % difference from
%   each curve's own finest N_cool_seg, so convergence is easy to read.
%
% For the station / HX profiles of ONE case across N_cool_seg (or models,
% N_segments, DNS), use compare_models.m.
%
% By default the input conditions are random ones drawn from the saved
% results (see design_source below).
%
% Every model evaluation is cached under results/id-<id>/ exactly like
% compare_models.m: if a matching results file exists it is loaded,
% otherwise the model is run and saved (overwrite = false, so nothing is
% ever clobbered and the first save of a case is not versioned).
% Re-running the script is therefore cheap.

clear; clc;

%% ---- Which input conditions to use ---- %%
% 'saved'  : pick n_cases random input conditions (ids) out of the ones
%            already stored in results/results_index.mat. Everything about
%            the condition (flight phase, temperatures, Q, e, r, d2, ...)
%            is taken from that id, so each case is exactly as it was run
%            before. Only ids from the discretised model (those that carry
%            tol_T) are eligible.
% 'manual' : the baseline of compare_models.m plus n_random random
%            (d2, e, r) draws, with every other input fixed (see below).
design_source = 'saved';
n_cases       = 4;        % 'saved': how many ids to draw
seed          = [];       % [] = different cases every run; a number = repeatable draw
n_random      = 3;        % 'manual': random draws on top of the baseline

%% ---- What to sweep ---- %%
N_cool_list = [1 2 3 5 7 10 15 20 30];   % must be >= 1 (0 would be the 1D path)
N_seg_list  = [50 20 10];                % first = finest; one line style each
line_styles = {'-', '--', ':', '-.'};    % matched to N_seg_list in order

use_DNS      = false;
save_results = true;
overwrite    = false;   % never clobber; the first save of a case is not versioned

%% ---- Build the list of input conditions (cases{i} = physical_inputs struct) ---- %%
if ~isempty(seed), rng(seed); else, rng('shuffle'); end

switch lower(design_source)
    case 'saved'
        idx_file = fullfile('results', 'results_index.mat');
        if ~isfile(idx_file)
            error('Ncoolseg_L_curves: %s not found -- run something with save_results = true first, or use design_source = ''manual''.', idx_file);
        end
        loaded = load(idx_file, 'results_index');
        ri = loaded.results_index;
        is_disc = cellfun(@(s) isfield(s, 'tol_T'), ri.inputs);
        elig = find(is_disc);
        if isempty(elig)
            error('Ncoolseg_L_curves: no discretised-model entries in the results index.');
        end
        pick = elig(randperm(numel(elig), min(n_cases, numel(elig))));
        cases   = ri.inputs(pick);
        case_id = ri.id(pick);
        fprintf('Picked %d of %d saved discretised input conditions: ids %s\n', ...
            numel(pick), numel(elig), mat2str(case_id(:)'));

    case 'manual'
        flight_phase = 3;
        base = physical_inputs(4 - 2*0.1/sqrt(2), 4 - 2*0.1/sqrt(2), 60, "OFF", 1, 44.4, 2, ...
            2.25e6, 70+273, 85+273, 4876, 22632, 216.65, 0.48, 3.8, 0.36, 128, 287, ...
            flight_phase, 0, 0, 0, 3);
        ranges.fr_dia    = [0.45 0.55];
        ranges.air_side  = [3 14];
        ranges.cool_side = [3 5];
        cases = cell(1 + n_random, 1);
        cases{1} = base;
        for i = 1:n_random
            c = base;
            c.d2_init = ranges.fr_dia(1)    + diff(ranges.fr_dia)   *rand;
            c.e       = ranges.air_side(1)  + diff(ranges.air_side) *rand;
            c.r       = ranges.cool_side(1) + diff(ranges.cool_side)*rand;
            cases{i+1} = c;
        end
        case_id = nan(numel(cases), 1);   % filled in once the cases are saved/looked up

    otherwise
        error('Ncoolseg_L_curves: design_source must be ''saved'' or ''manual''.');
end
n_dp = numel(cases);

%% ---- Run / load every combination ---- %%
n_Nc = numel(N_cool_list);
n_Ns = numel(N_seg_list);
L = nan(n_dp, n_Ns, n_Nc);

for i = 1:n_dp
    S = cases{i};
    for j = 1:n_Ns
        for k = 1:n_Nc
            fprintf('design %d/%d  N_seg=%d  N_cool=%d ... ', i, n_dp, N_seg_list(j), N_cool_list(k));
            try
                res = get_or_run(use_DNS, N_seg_list(j), N_cool_list(k), S, save_results, overwrite);
                L(i,j,k) = res.L_solution;
                fprintf('L = %.4f m\n', res.L_solution);
            catch ME
                fprintf('FAILED: %s\n', ME.message);
            end
        end
    end
end

cols = lines(n_dp);

figure('Name', 'L_solution vs N_cool_seg', 'NumberTitle', 'off'); hold on; grid on
for i = 1:n_dp
    for j = 1:n_Ns
        plot(N_cool_list, squeeze(L(i,j,:)), line_styles{mod(j-1,numel(line_styles))+1}, ...
            'Color', cols(i,:), 'LineWidth', 1.4, 'Marker', 'o', 'MarkerSize', 4, ...
            'MarkerFaceColor', cols(i,:), 'HandleVisibility', 'off');
    end
end
% Legend proxies: colours = design points, line styles = N_segments
for i = 1:n_dp
    plot(nan, nan, '-', 'Color', cols(i,:), 'LineWidth', 2, ...
        'DisplayName', case_label(cases{i}, case_id(i)));
end
for j = 1:n_Ns
    plot(nan, nan, line_styles{mod(j-1,numel(line_styles))+1}, 'Color', [0.3 0.3 0.3], ...
        'LineWidth', 1.5, 'DisplayName', sprintf('N_{segments} = %d', N_seg_list(j)));
end
xlabel('N_{cool,seg}'); ylabel('L_{solution} [m]');
title('HX length vs. coolant-side resolution'); legend('Location', 'best');

% Same curves relative to each curve's own finest N_cool_seg
figure('Name', 'L_solution vs N_cool_seg (relative)', 'NumberTitle', 'off'); hold on; grid on
for i = 1:n_dp
    for j = 1:n_Ns
        y = squeeze(L(i,j,:));
        yref = y(find(~isnan(y), 1, 'last'));
        plot(N_cool_list, 100*(y - yref)/yref, line_styles{mod(j-1,numel(line_styles))+1}, ...
            'Color', cols(i,:), 'LineWidth', 1.4, 'Marker', 'o', 'MarkerSize', 4, ...
            'MarkerFaceColor', cols(i,:), 'HandleVisibility', 'off');
    end
end
yline(0, 'k-', 'HandleVisibility', 'off');
yline([-1 1], 'k:', 'HandleVisibility', 'off');
xlabel('N_{cool,seg}'); ylabel(sprintf('L_{solution} vs. value at N_{cool,seg} = %d [%%]', N_cool_list(end)));
title('Convergence of HX length (dotted lines = \pm1 %)');

%% ---- Local helper: legend text for one input condition ---- %%
function s = case_label(S, id)
if isnan(id)
    id = find_matching_id(S);
end
if isempty(id) || isnan(id), idtxt = 'new'; else, idtxt = sprintf('id-%d', id); end
s = sprintf('%s: d_2=%.3f, e=%.2f, r=%.2f, phase %d', idtxt, S.d2_init, S.e, S.r, S.flight_phase);
end

%% ---- Local helper: load from cache, else run + save ---- %%
% S is a physical_inputs struct; it is itself the results_index key, so the
% lookup matches exactly the id the stored results were saved under.
function res = get_or_run(use_DNS, N_seg, N_cool, S, save_results, overwrite)
id = find_matching_id(S);

res = [];
if ~isempty(id)
    f = fullfile('results', sprintf('id-%d', id), ...
        disc_results_filename(N_seg, N_cool, use_DNS, id, S.solve_T));
    if exist(f, 'file')
        loaded = load(f, 'results');
        res = loaded.results;
    end
end

if isempty(res)
    res = run_disc_model_fwdpass(use_DNS, N_seg, N_cool, S.e, S.r, S.hx_theta, S.fan, S.fpr_init, ...
        S.M_dot_coolant, S.n_modules, S.Q_tot, S.T_in_fc, S.T_out_fc, S.h, S.p11, S.t11, ...
        S.d2_init, S.AR_diff, S.AR_noz, S.V_inf, S.R, S.flight_phase, S.M_dot_FOD, S.M_dot_comp, ...
        S.solve_T, save_results, overwrite, S.tol_T);
end
end