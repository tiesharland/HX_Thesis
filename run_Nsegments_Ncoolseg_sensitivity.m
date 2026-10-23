%% Nsegments_Ncoolseg_sensitivity_study.m
% Monte Carlo convergence study of
%   Part A: N_segments (air-side resolution), 1D path (N_cool_seg = 0)
%   Part B: N_cool_seg (coolant-side resolution, 2D path), repeated at a
%           few N_segments values
% for every model variant in use_DNS_list (correlation-based and/or DNS).
%
% SAVING (overnight-safe):
%   * Every model evaluation is recomputed (nothing is read from cache) and
%     saved to results/id-<id>/ with overwrite = false, i.e. a rerun adds a
%     new version (a50_<id>-1.mat, -2, ...) and never clobbers an old file.
%   * The samples are processed in chunks. After each chunk the evaluations
%     are saved serially on the client (the parfor workers never touch
%     results_index.mat) and then dropped from memory, so a crash only
%     loses the current chunk and memory stays flat.
%   * The study summary .mat files get a timestamp in their name, so reruns
%     never overwrite them either.
%   * DNS: the thermoturb query log is reset before every model run, so each
%     saved result only carries its own queries (otherwise the log would
%     keep growing on every worker and bloat each saved file).

clear; clc;

%% ---- Fixed inputs (mirroring Copy_of_tmsObjective) ---- %%
AR_noz          = 0.36;
AR_diff         = 3.8;
n_modules       = 2;
M_dot_coolant   = 44.4;
hx_theta        = 60;
fpr_init        = 1;
fan             = "OFF";
M_dot_FOD       = 0;
M_dot_comp      = 0;
tol_T           = 0;
solve_T         = 0;        % bookkeeping only, not yet implemented -- always 0

use_DNS_list    = [false true];   % model variants to run, in this order (same design points for each)

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

p11 = 22632;
t11 = 216.65;
R   = 287;

% Everything run_disc_model_fwdpass / save_disc_results need besides the
% per-evaluation (use_DNS, N_segments, N_cool_seg, e, r, d2).
P = struct('hx_theta',hx_theta, 'fan',fan, 'fpr_init',fpr_init, ...
    'M_dot_coolant',M_dot_coolant, 'n_modules',n_modules, 'Q_tot',Q_tot, ...
    'T_in_fc',T_in_fc, 'T_out_fc',T_out_fc, 'h',h, 'p11',p11, 't11',t11, ...
    'AR_diff',AR_diff, 'AR_noz',AR_noz, 'V_inf',V_inf, 'R',R, ...
    'flight_phase',flight_phase, 'M_dot_FOD',M_dot_FOD, 'M_dot_comp',M_dot_comp, ...
    'solve_T',solve_T, 'tol_T',tol_T);

%% ---- Design variable ranges ---- %%
ranges.fr_dia    = [0.45 0.55];
ranges.air_side  = [3 8];
ranges.cool_side = [3 5];

%% ---- Candidates + convergence tolerance ---- %%
N_segments_candidates = [1 5 10 25 50 100];
N_cool_seg_fixed      = 0;

N_cool_seg_candidates = [1 2 5 10 15 20 30];
N_segments_for_coolseg_study = [5 10 25 50];

tol          = 0.01;  % 1%
metric_names = {'drag_tot', 'L_solution', 'dp_hx'};

%% ---- Monte Carlo sampling ---- %%
n_samples  = 500;
n_workers  = 10;                 % parallel workers (laptop has 16 processors: leave a few free,
                                 % and each DNS worker also holds its own copy of the table)
chunk_size = 3*n_workers;        % samples per parfor batch; saved + freed after each batch

p = gcp('nocreate');
if isempty(p)
    p = parpool(n_workers);
elseif p.NumWorkers ~= n_workers
    fprintf('Note: an existing pool with %d workers is open (n_workers = %d requested).\n', p.NumWorkers, n_workers);
end
% The DNS table (thermoturb_table.mat) and the model files are found
% relative to the current folder, so make sure every worker is in it.
pctRunOnAll(sprintf('cd(''%s'')', strrep(pwd, '''', '''''')));

rng(1);
design_points = table( ...
    ranges.fr_dia(1)    + diff(ranges.fr_dia)   *rand(n_samples,1), ...
    ranges.air_side(1)  + diff(ranges.air_side) *rand(n_samples,1), ...
    ranges.cool_side(1) + diff(ranges.cool_side)*rand(n_samples,1), ...
    'VariableNames', {'fr_dia','air_side','cool_side'});

stamp = datestr(now, 'yyyymmdd_HHMMSS');
t_start = tic;

for v = 1:numel(use_DNS_list)
use_DNS = use_DNS_list(v);
if use_DNS, tag = 'dns'; else, tag = 'corr'; end
fprintf('\n################ Variant: %s (use_DNS = %d) ################\n', tag, use_DNS);

%% ==================================================================== %%
%% ---- PART A: N_segments (air-side resolution) convergence study ---- %%
%% ==================================================================== %%

N_conv_all      = nan(n_samples,1);
N_conv_drag     = nan(n_samples,1);
N_conv_length   = nan(n_samples,1);
N_conv_dp       = nan(n_samples,1);
sample_status   = strings(n_samples,1); % "ok" | "no_convergence" | "infeasible_at_N" | "sample_error"
n_saved = 0;

for c0 = 1:chunk_size:n_samples
    idx = c0:min(c0+chunk_size-1, n_samples);
    n_c = numel(idx);
    c_all = nan(n_c,1); c_drag = nan(n_c,1); c_len = nan(n_c,1); c_dp = nan(n_c,1);
    c_status = strings(n_c,1);
    c_entries = cell(n_c,1);

    parfor ii = 1:n_c
        i  = idx(ii);
        d2 = design_points.fr_dia(i);
        e  = design_points.air_side(i);
        r  = design_points.cool_side(i);

        collector = ResultsCollector();
        model_fun = @(N) capture(run_model(use_DNS, N, N_cool_seg_fixed, e, r, d2, P), ...
            collector, e, r, d2, N, N_cool_seg_fixed);

        try
            [N_conv, N_conv_per_metric, ~, status] = find_converged_Nsegments(model_fun, N_segments_candidates, tol, metric_names);
            c_all(ii)    = N_conv;
            c_drag(ii)   = N_conv_per_metric.drag_tot;
            c_len(ii)    = N_conv_per_metric.L_solution;
            c_dp(ii)     = N_conv_per_metric.dp_hx;
            c_status(ii) = status;
        catch ME
            warning('Sample %d failed entirely (outside find_converged_Nsegments): %s', i, ME.message);
            c_status(ii) = "sample_error";
        end

        c_entries{ii} = collector.entries;
    end

    N_conv_all(idx) = c_all;  N_conv_drag(idx) = c_drag;
    N_conv_length(idx) = c_len;  N_conv_dp(idx) = c_dp;
    sample_status(idx) = c_status;

    % Serial save of this chunk (overwrite = false: new version, never clobbers)
    n_saved = n_saved + save_entries(c_entries, P, use_DNS);
    clear c_entries
    fprintf('[%s] Part A: %d/%d samples done, %d evaluations saved (%.1f min elapsed)\n', ...
        tag, idx(end), n_samples, n_saved, toc(t_start)/60);
end

%% ---- Failure reporting (before filtering to 'valid') ---- %%
n_ok             = sum(sample_status == "ok");
n_no_convergence = sum(sample_status == "no_convergence");
n_infeasible     = sum(sample_status == "infeasible_at_N");
n_sample_error   = sum(sample_status == "sample_error");

fprintf('\n=== [%s] Part A: N_segments sample outcome breakdown (%d total design points) ===\n', tag, n_samples);
fprintf('  Converged within N_candidates      : %d (%.0f%%)\n', n_ok, 100*n_ok/n_samples);
fprintf('  Ran, full sweep, never converged   : %d (%.0f%%)\n', n_no_convergence, 100*n_no_convergence/n_samples);
fprintf('  Infeasible design (failed at some N): %d (%.0f%%)\n', n_infeasible, 100*n_infeasible/n_samples);
fprintf('  Failed entirely (unexpected error) : %d (%.0f%%)\n', n_sample_error, 100*n_sample_error/n_samples);

n_excluded = n_no_convergence + n_infeasible + n_sample_error;
if n_excluded > 0
    fprintf(['  NOTE: %.0f%% of sampled design points are excluded from the ' ...
             'statistics below (discarded as infeasible or otherwise unusable). ' ...
             'The reported distribution describes only the feasible/converged subset.\n'], ...
             100*n_excluded/n_samples);
end

valid = (sample_status == "ok");

% Saved right away (before plotting) so a plotting error can't lose it.
save(sprintf('Nsegments_convergence_study_%s_%s.mat', tag, stamp), 'design_points', 'N_conv_all', ...
     'N_conv_drag', 'N_conv_length', 'N_conv_dp', 'N_segments_candidates', 'tol', 'sample_status', 'use_DNS');

if any(valid)
    figure('Name',sprintf('[%s] Part A: Feasibility across design space',tag),'NumberTitle','off');
    gscatter(design_points.fr_dia, design_points.air_side, sample_status);
    xlabel('fr\_dia'); ylabel('air\_side');
    title(sprintf('[%s] Part A (N_{segments} sweep): sample outcome across sampled design space', tag));

    fprintf('\n=== [%s] N_segments convergence across %d FEASIBLE design points ===\n', tag, sum(valid));
    fprintf('Overall (all metrics simultaneous, tol=%.1f%%):\n', tol*100);
    fprintf('  Min    = %.1f\n', min(N_conv_all(valid)));
    fprintf('  Mean   = %.1f\n', mean(N_conv_all(valid)));
    fprintf('  Median = %.1f\n', median(N_conv_all(valid)));
    fprintf('  Max    = %.1f\n', max(N_conv_all(valid)));
    fprintf('  Std    = %.1f\n', std(N_conv_all(valid)));

    fprintf('\nPer-metric medians (diagnostic -- which metric is slowest to converge):\n');
    fprintf('  drag_tot    : median N = %.1f  (range %.0f-%.0f)\n', ...
            median(N_conv_drag(valid),'omitnan'), min(N_conv_drag(valid)), max(N_conv_drag(valid)));
    fprintf('  L_solution  : median N = %.1f  (range %.0f-%.0f)\n', ...
            median(N_conv_length(valid),'omitnan'), min(N_conv_length(valid)), max(N_conv_length(valid)));
    fprintf('  dp_hx       : median N = %.1f  (range %.0f-%.0f)\n', ...
            median(N_conv_dp(valid),'omitnan'), min(N_conv_dp(valid)), max(N_conv_dp(valid)));

    figure('Name',sprintf('[%s] Part A: N_segments Convergence Distribution',tag),'NumberTitle','off');
    subplot(1,2,1)
    histogram(N_conv_all(valid), 'BinMethod','integers', 'FaceColor',[0.839 0.153 0.157]);
    xlabel('Converged N_{segments}'); ylabel('Count');
    title(sprintf('[%s] Distribution across %d feasible design points', tag, sum(valid)));
    grid on
    subplot(1,2,2)
    boxplot([N_conv_drag(valid), N_conv_length(valid), N_conv_dp(valid)], ...
            'Labels', {'drag\_tot','L\_solution','dp\_hx'});
    ylabel('Converged N_{segments}');
    title('Per-metric convergence N (which metric is slowest?)');
    grid on
else
    warning('[%s] Part A: no converged design points -- skipping statistics/plots.', tag);
end

%% ==================================================================== %%
%% ---- PART B: N_cool_seg (coolant-side resolution) convergence   ---- %%
%% ----         study, repeated across a few N_segments values     ---- %%
%% ==================================================================== %%

n_outer = numel(N_segments_for_coolseg_study);

N_conv_coolseg_all      = nan(n_samples, n_outer);
N_conv_coolseg_drag     = nan(n_samples, n_outer);
N_conv_coolseg_length   = nan(n_samples, n_outer);
N_conv_coolseg_dp       = nan(n_samples, n_outer);
sample_status_coolseg   = strings(n_samples, n_outer);

for j = 1:n_outer
    N_seg_fixed = N_segments_for_coolseg_study(j);
    fprintf('\n--- [%s] Part B: N_cool_seg convergence study at N_segments = %d (%d/%d) ---\n', ...
            tag, N_seg_fixed, j, n_outer);
    n_saved_j = 0;

    for c0 = 1:chunk_size:n_samples
        idx = c0:min(c0+chunk_size-1, n_samples);
        n_c = numel(idx);
        c_all = nan(n_c,1); c_drag = nan(n_c,1); c_len = nan(n_c,1); c_dp = nan(n_c,1);
        c_status = strings(n_c,1);
        c_entries = cell(n_c,1);

        parfor ii = 1:n_c
            i  = idx(ii);
            d2 = design_points.fr_dia(i);
            e  = design_points.air_side(i);
            r  = design_points.cool_side(i);

            collector = ResultsCollector();
            model_fun = @(Nc) capture(run_model(use_DNS, N_seg_fixed, Nc, e, r, d2, P), ...
                collector, e, r, d2, N_seg_fixed, Nc);

            try
                [N_conv, N_conv_per_metric, ~, status] = find_converged_Nsegments(model_fun, N_cool_seg_candidates, tol, metric_names);
                c_all(ii)    = N_conv;
                c_drag(ii)   = N_conv_per_metric.drag_tot;
                c_len(ii)    = N_conv_per_metric.L_solution;
                c_dp(ii)     = N_conv_per_metric.dp_hx;
                c_status(ii) = status;
            catch ME
                warning('Sample %d (N_segments=%d) failed entirely: %s', i, N_seg_fixed, ME.message);
                c_status(ii) = "sample_error";
            end

            c_entries{ii} = collector.entries;
        end

        N_conv_coolseg_all(idx,j) = c_all;  N_conv_coolseg_drag(idx,j) = c_drag;
        N_conv_coolseg_length(idx,j) = c_len;  N_conv_coolseg_dp(idx,j) = c_dp;
        sample_status_coolseg(idx,j) = c_status;

        n_saved_j = n_saved_j + save_entries(c_entries, P, use_DNS);
        clear c_entries
        fprintf('[%s] Part B (N_seg=%d): %d/%d samples done, %d evaluations saved (%.1f min elapsed)\n', ...
            tag, N_seg_fixed, idx(end), n_samples, n_saved_j, toc(t_start)/60);
    end

    % Checkpoint after each N_segments value (timestamped, never overwritten)
    save(sprintf('Ncoolseg_convergence_study_%s_%s.mat', tag, stamp), 'design_points', 'N_conv_coolseg_all', ...
         'N_conv_coolseg_drag', 'N_conv_coolseg_length', 'N_conv_coolseg_dp', ...
         'N_cool_seg_candidates', 'N_segments_for_coolseg_study', 'tol', 'sample_status_coolseg', 'use_DNS');

    status_j = sample_status_coolseg(:,j);
    n_ok_j             = sum(status_j == "ok");
    n_no_convergence_j = sum(status_j == "no_convergence");
    n_infeasible_j     = sum(status_j == "infeasible_at_N");
    n_sample_error_j   = sum(status_j == "sample_error");

    fprintf('  Converged within N_cool_seg_candidates : %d (%.0f%%)\n', n_ok_j, 100*n_ok_j/n_samples);
    fprintf('  Ran, full sweep, never converged       : %d (%.0f%%)\n', n_no_convergence_j, 100*n_no_convergence_j/n_samples);
    fprintf('  Infeasible design (failed at some N)    : %d (%.0f%%)\n', n_infeasible_j, 100*n_infeasible_j/n_samples);
    fprintf('  Failed entirely (unexpected error)      : %d (%.0f%%)\n', n_sample_error_j, 100*n_sample_error_j/n_samples);

    valid_j = (status_j == "ok");

    figure('Name', sprintf('[%s] Part B: Feasibility (N_segments=%d)', tag, N_seg_fixed), 'NumberTitle', 'off');
    gscatter(design_points.fr_dia, design_points.cool_side, status_j);
    xlabel('fr\_dia'); ylabel('cool\_side');
    title(sprintf('[%s] Part B (N_{cool,seg} sweep, N_{segments}=%d): sample outcome', tag, N_seg_fixed));

    if any(valid_j)
        fprintf('  Converged N_cool_seg (feasible subset, n=%d): median = %.1f, range %.0f-%.0f\n', ...
                sum(valid_j), median(N_conv_coolseg_all(valid_j,j)), ...
                min(N_conv_coolseg_all(valid_j,j)), max(N_conv_coolseg_all(valid_j,j)));

        figure('Name', sprintf('[%s] Part B: N_cool_seg Convergence Distribution (N_segments=%d)', tag, N_seg_fixed), 'NumberTitle', 'off');
        subplot(1,2,1)
        histogram(N_conv_coolseg_all(valid_j,j), 'BinMethod','integers', 'FaceColor',[0.122 0.467 0.706]);
        xlabel('Converged N_{cool,seg}'); ylabel('Count');
        title(sprintf('[%s] N_{segments}=%d: %d feasible points', tag, N_seg_fixed, sum(valid_j)));
        grid on
        subplot(1,2,2)
        boxplot([N_conv_coolseg_drag(valid_j,j), N_conv_coolseg_length(valid_j,j), N_conv_coolseg_dp(valid_j,j)], ...
                'Labels', {'drag\_tot','L\_solution','dp\_hx'});
        ylabel('Converged N_{cool,seg}');
        title('Per-metric convergence N (which metric is slowest?)');
        grid on
    else
        warning('[%s] Part B (N_segments=%d): no feasible/converged design points -- skipping distribution plot.', tag, N_seg_fixed);
    end
end

%% ---- Cross-comparison: does converged N_cool_seg depend on N_segments? ---- %%
if any(any(sample_status_coolseg == "ok", 1))
    figure('Name',sprintf('[%s] Part B: N_cool_seg convergence vs. air-side resolution',tag),'NumberTitle','off');
    x = []; g = [];
    for j = 1:n_outer
        mask = sample_status_coolseg(:,j) == "ok";
        x = [x; N_conv_coolseg_all(mask, j)]; %#ok<AGROW>
        g = [g; repmat(N_segments_for_coolseg_study(j), nnz(mask), 1)]; %#ok<AGROW>
    end
    boxplot(x, g);
    xlabel('N_{segments} (air-side resolution)'); ylabel('Converged N_{cool,seg}');
    title(sprintf('[%s] Does the coolant-side resolution needed depend on the air-side resolution?', tag));
    grid on
end

drawnow;
end   % variants

fprintf('\nDone. Total time %.1f h.\n', toc(t_start)/3600);

%% ---- Local helper: one model evaluation (no saving; the client saves) ---- %%
% For DNS, the thermoturb query log is reset first so this result's saved
% log holds only its own queries.
function res = run_model(use_DNS, N_segments, N_cool_seg, e, r, d2, P)
if use_DNS
    thermoturb_query_log('reset');
end
res = run_disc_model_fwdpass(use_DNS, N_segments, N_cool_seg, e, r, P.hx_theta, P.fan, P.fpr_init, ...
    P.M_dot_coolant, P.n_modules, P.Q_tot, P.T_in_fc, P.T_out_fc, P.h, P.p11, P.t11, ...
    d2, P.AR_diff, P.AR_noz, P.V_inf, P.R, P.flight_phase, P.M_dot_FOD, P.M_dot_comp, ...
    P.solve_T, false, false, P.tol_T);
end

%% ---- Local helper: tag a model evaluation into this iteration's collector ---- %%
% Returned untouched, so find_converged_Nsegments sees exactly what
% run_model returned; the side effect (collector.add) is what makes every
% evaluation recoverable for the serial save afterwards (ResultsCollector.m).
function res = capture(res, collector, e, r, d2, N_segments, N_cool_seg)
collector.add(res, e, r, d2, N_segments, N_cool_seg);
end

%% ---- Local helper: save a chunk's evaluations serially ---- %%
% overwrite = false -> a rerun saves new versions, never clobbers.
function n = save_entries(entries_cell, P, use_DNS)
n = 0;
for i = 1:numel(entries_cell)
    for k = 1:numel(entries_cell{i})
        ent = entries_cell{i}{k};
        save_disc_results(ent.results, ent.e, ent.r, P.hx_theta, P.fan, P.fpr_init, ...
            P.M_dot_coolant, P.n_modules, P.Q_tot, P.T_in_fc, P.T_out_fc, P.h, P.p11, P.t11, ...
            ent.d2, P.AR_diff, P.AR_noz, P.V_inf, P.R, P.flight_phase, P.M_dot_FOD, P.M_dot_comp, ...
            P.solve_T, P.tol_T, ent.N_segments, ent.N_cool_seg, use_DNS, false);
        n = n + 1;
    end
end
end