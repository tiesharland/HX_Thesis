%% run_Ncoolseg_L_sensitivity.m
% How does the sized HX length L_solution change with N_cool_seg (the
% coolant-side resolution of the 2D direct-marching path)?
%
%   One COLOUR per input condition (id), one LINE STYLE per N_segments
%   (solid = finest). A second figure shows the same curves as % difference
%   from each curve's own finest N_cool_seg, so convergence is easy to read.
%
% This script NEVER runs the model. It only reads results that are already
% saved under results/id-<id>/ (the files made by run_disc_model_fwdpass /
% the sensitivity study / compare_models). Combinations that were never run
% are simply left out of the curves.
%
% For the station / HX profiles of ONE case across N_cool_seg (or models,
% N_segments, DNS), use compare_models.m.

clear; clc;

%% ---- Settings ---- %%
ids         = [];       % specific ids to plot, e.g. [3 7 12]; [] = draw randomly (below)
n_cases     = 4;        % random draw: how many saved ids
seed        = [];        % random draw: [] = different ids every run; a number = repeatable
use_DNS     = false;    % which saved variant to read (correlation-based vs DNS-based)

N_cool_list = [1 2 3 5 7 10 15 20 30];   % must be >= 1 (0 would be the 1D path)
N_seg_list  = [50 20 10];                % first = finest; one line style each
line_styles = {'-', '--', ':', '-.'};    % matched to N_seg_list in order

min_points  = 2;        % random draw only: an id needs at least this many saved N_cool_seg
                        % values at some N_segments to be eligible (one point is not a curve)

%% ---- Find the saved ids that have usable results ---- %%
idx_file = fullfile('results', 'results_index.mat');
if ~isfile(idx_file)
    error('run_Ncoolseg_L_sensitivity: %s not found -- run the model with save_results = true first.', idx_file);
end
loaded = load(idx_file, 'results_index');
ri = loaded.results_index;

n_Nc = numel(N_cool_list);
n_Ns = numel(N_seg_list);

is_disc = cellfun(@(s) isfield(s, 'tol_T'), ri.inputs);   % discretised-model ids only

% saved_file{row}{j,k} = path of the saved result for N_seg_list(j), N_cool_list(k), or ''
saved_file = cell(height(ri), 1);
for row = 1:height(ri)
    if ~is_disc(row), continue, end
    if ~isempty(ids) && ~ismember(ri.id(row), ids), continue, end   % only scan what is needed
    id = ri.id(row);
    S  = ri.inputs{row};
    files = repmat({''}, n_Ns, n_Nc);
    for j = 1:n_Ns
        for k = 1:n_Nc
            f = fullfile('results', sprintf('id-%d', id), ...
                disc_results_filename(N_seg_list(j), N_cool_list(k), use_DNS, id, S.solve_T));
            if isfile(f), files{j,k} = f; end
        end
    end
    saved_file{row} = files;
end

if ~isempty(ids)
    % Chosen ids, in the order given
    pick = zeros(0,1);
    for id = ids(:)'
        row = find(ri.id == id & is_disc, 1);
        if isempty(row)
            warning('id %d is not a discretised-model id in the results index -- skipped.', id);
        elseif ~any(cellfun(@(x) ~isempty(x), saved_file{row}), 'all')
            warning('id-%d has no saved results for these N_cool_list/N_seg_list/use_DNS -- skipped.', id);
        else
            pick(end+1,1) = row; %#ok<SAGROW>
        end
    end
    if isempty(pick)
        error('run_Ncoolseg_L_sensitivity: none of the requested ids has usable saved results.');
    end
    fprintf('Using requested ids: %s\n', mat2str(ri.id(pick)'));
else
    % Random draw among ids with at least min_points saved N_cool_seg values at some N_segments
    elig = [];
    for row = find(is_disc)'
        have = cellfun(@(x) ~isempty(x), saved_file{row});
        if any(sum(have, 2) >= min_points), elig(end+1) = row; end %#ok<SAGROW>
    end
    if isempty(elig)
        error(['run_Ncoolseg_L_sensitivity: no saved id has >= %d N_cool_seg results at the same ' ...
               'N_segments (N_cool_list/N_seg_list/use_DNS as set above).'], min_points);
    end
    if ~isempty(seed), rng(seed); else, rng('shuffle'); end
    pick = elig(randperm(numel(elig), min(n_cases, numel(elig))));
    fprintf('Picked %d of %d eligible saved ids: %s\n', numel(pick), numel(elig), mat2str(ri.id(pick)'));
end

cases   = ri.inputs(pick);
case_id = ri.id(pick);
n_dp    = numel(pick);

%% ---- Load the saved results ---- %%
L = nan(n_dp, n_Ns, n_Nc);
for i = 1:n_dp
    files = saved_file{pick(i)};
    for j = 1:n_Ns
        for k = 1:n_Nc
            if ~isempty(files{j,k})
                loaded = load(files{j,k}, 'results');
                L(i,j,k) = loaded.results.L_solution;
            end
        end
    end
    fprintf('id-%d: %d of %d combinations saved\n', case_id(i), nnz(~isnan(L(i,:,:))), n_Ns*n_Nc);
end

%% ---- Plots ---- %%
cols = lines(n_dp);

figure('Name', 'L_solution vs N_cool_seg', 'NumberTitle', 'off'); hold on; grid on
for i = 1:n_dp
    for j = 1:n_Ns
        y  = squeeze(L(i,j,:))';
        ok = ~isnan(y);
        if ~any(ok), continue, end
        plot(N_cool_list(ok), y(ok), line_styles{mod(j-1,numel(line_styles))+1}, ...
            'Color', cols(i,:), 'LineWidth', 1.4, 'Marker', 'o', 'MarkerSize', 4, ...
            'MarkerFaceColor', cols(i,:), 'HandleVisibility', 'off');
    end
end
% Legend proxies: colours = ids, line styles = N_segments
for i = 1:n_dp
    plot(nan, nan, '-', 'Color', cols(i,:), 'LineWidth', 2, 'DisplayName', case_label(cases{i}, case_id(i)));
end
for j = 1:n_Ns
    plot(nan, nan, line_styles{mod(j-1,numel(line_styles))+1}, 'Color', [0.3 0.3 0.3], ...
        'LineWidth', 1.5, 'DisplayName', sprintf('N_{segments} = %d', N_seg_list(j)));
end
xlabel('N_{cool,seg}'); ylabel('L_{solution} [m]');
title('HX length vs. coolant-side resolution (saved results)'); legend('Location', 'best');

% Same curves relative to each curve's own finest saved N_cool_seg
figure('Name', 'L_solution vs N_cool_seg (relative)', 'NumberTitle', 'off'); hold on; grid on
for i = 1:n_dp
    for j = 1:n_Ns
        y  = squeeze(L(i,j,:))';
        ok = ~isnan(y);
        if nnz(ok) < 1, continue, end
        yref = y(find(ok, 1, 'last'));
        plot(N_cool_list(ok), 100*(y(ok) - yref)/yref, line_styles{mod(j-1,numel(line_styles))+1}, ...
            'Color', cols(i,:), 'LineWidth', 1.4, 'Marker', 'o', 'MarkerSize', 4, ...
            'MarkerFaceColor', cols(i,:), 'HandleVisibility', 'off');
    end
end
yline(0, 'k-', 'HandleVisibility', 'off');
yline([-1 1], 'k:', 'HandleVisibility', 'off');
xlabel('N_{cool,seg}'); ylabel('L_{solution} vs. finest saved N_{cool,seg} [%]');
title('Convergence of HX length (dotted lines = \pm1 %)');

%% ---- Local helper: legend text for one input condition ---- %%
function s = case_label(S, id)
s = sprintf('id-%d: d_2=%.3f, e=%.2f, r=%.2f, phase %d', id, S.d2_init, S.e, S.r, S.flight_phase);
end