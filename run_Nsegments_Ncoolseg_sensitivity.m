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
use_DNS         = false;    % correlation-based: fast enough for a broad sweep like this

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

%% ---- Design variable ranges ---- %%
ranges.fr_dia    = [0.45 0.55];
ranges.air_side  = [3 14];
ranges.cool_side = [3 5];

%% ---- Candidates + convergence tolerance ---- %%
% Part A (N_segments, air-side resolution): 1D iterative path, i.e.
% N_cool_seg fixed at 0 (HX_design1_disc's "old" variant) -- matches how
% this study was always run before N_cool_seg/the 2D path existed.
N_segments_candidates = [1 2 5 10 15 20 30 40 50 75 100 150 200 250];
N_cool_seg_fixed       = 0;

% Part B (N_cool_seg, coolant-side resolution of the 2D direct-marching
% path): N_cool_seg must be >= 1 to even select that path (0 means the 1D
% path instead, which has no coolant-side discretisation to be sensitive
% to) -- see HX_design1_disc.m / run_disc_model_fwdpass.m's docstrings.
N_cool_seg_candidates = [1 2 3 5 7 10 15 20 30];

% Repeat Part B across a few representative air-side resolutions, to see
% whether the converged N_cool_seg depends on how finely the air side
% happens to be discretised. Kept short (not a full sweep in its own
% right) since Part B's cost scales with n_samples * numel(this) *
% numel(N_cool_seg_candidates) -- trim this list (or N_cool_seg_candidates,
% or n_samples below) first if runtime becomes a problem.
N_segments_for_coolseg_study = [10 50 100];

tol          = 0.01;  % 1%
metric_names = {'drag_tot', 'L_solution', 'dp_hx'};

%% ---- Monte Carlo sampling ---- %%
% Shared between Part A and Part B so a given design point's air-side and
% coolant-side convergence behaviour can be compared directly.
n_samples = 500;

p = gcp('nocreate');
if isempty(p)
    parpool(4);
end

rng(1);
design_points = table( ...
    ranges.fr_dia(1)    + diff(ranges.fr_dia)   *rand(n_samples,1), ...
    ranges.air_side(1)  + diff(ranges.air_side) *rand(n_samples,1), ...
    ranges.cool_side(1) + diff(ranges.cool_side)*rand(n_samples,1), ...
    'VariableNames', {'fr_dia','air_side','cool_side'});

%% ==================================================================== %%
%% ---- PART A: N_segments (air-side resolution) convergence study ---- %%
%% ==================================================================== %%

N_conv_all      = nan(n_samples,1);
N_conv_drag     = nan(n_samples,1);
N_conv_length   = nan(n_samples,1);
N_conv_dp       = nan(n_samples,1);
sample_status   = strings(n_samples,1); % "ok" | "no_convergence" | "infeasible_at_N" | "sample_error"
all_entries     = cell(n_samples,1);    % every model evaluation made, for the post-loop save pass

parfor i = 1:n_samples
    d2 = design_points.fr_dia(i);
    e  = design_points.air_side(i);
    r  = design_points.cool_side(i);

    % Every model evaluation this makes (one per N in N_segments_candidates,
    % up to the point of failure/infeasibility) gets stashed in collector
    % -- see ResultsCollector.m for why this (rather than save_results=true
    % here) is what keeps this parfor loop safe: saving to the SHARED
    % results/results_index.mat happens in a serial pass below instead,
    % after every worker's compute is done, so nobody races to update it.
    collector = ResultsCollector();
    model_fun = @(N) capture( ...
        run_disc_model_fwdpass(use_DNS, N, N_cool_seg_fixed, e, r, hx_theta, fan, fpr_init, ...
            M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
            d2, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, ...
            solve_T, false, false, tol_T), ...
        collector, e, r, d2, N, N_cool_seg_fixed);

    try
        [N_conv, N_conv_per_metric, ~, status] = find_converged_Nsegments(model_fun, N_segments_candidates, tol, metric_names);
        N_conv_all(i)    = N_conv;
        N_conv_drag(i)   = N_conv_per_metric.drag_tot;
        N_conv_length(i) = N_conv_per_metric.L_solution;
        N_conv_dp(i)     = N_conv_per_metric.dp_hx;
        sample_status(i) = status;
    catch ME
        warning('Sample %d failed entirely (outside find_converged_Nsegments): %s', i, ME.message);
        sample_status(i) = "sample_error";
    end

    all_entries{i} = collector.entries;
end

%% ---- Save every Part A model evaluation (serially -- see ResultsCollector.m) ---- %%
fprintf('\nSaving Part A (N_segments sweep) model evaluations to disk...\n');
n_saved = 0;
for i = 1:n_samples
    for k = 1:numel(all_entries{i})
        ent = all_entries{i}{k};
        save_disc_results(ent.results, ent.e, ent.r, hx_theta, fan, fpr_init, ...
            M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
            ent.d2, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, ...
            solve_T, tol_T, ent.N_segments, ent.N_cool_seg, use_DNS, false);   % overwrite=false: new version, never clobbers
        n_saved = n_saved + 1;
    end
end
fprintf('Saved %d model evaluations under results/ (new versions throughout, nothing overwritten).\n', n_saved);

%% ---- Failure reporting (before filtering to 'valid') ---- %%
n_ok             = sum(sample_status == "ok");
n_no_convergence = sum(sample_status == "no_convergence");
n_infeasible     = sum(sample_status == "infeasible_at_N");
n_sample_error   = sum(sample_status == "sample_error");

fprintf('\n=== Part A: N_segments sample outcome breakdown (%d total design points) ===\n', n_samples);
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

%% ---- Feasibility across design space (diagnostic) ---- %%
figure('Name','Part A: Feasibility across design space','NumberTitle','off');
gscatter(design_points.fr_dia, design_points.air_side, sample_status);
xlabel('fr\_dia'); ylabel('air\_side');
title('Part A (N_{segments} sweep): sample outcome across sampled design space');

%% ---- Summary statistics (feasible subset only) ---- %%
fprintf('\n=== N_segments convergence across %d FEASIBLE design points ===\n', sum(valid));
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

%% ---- Plots ---- %%
figure('Name','Part A: N_segments Convergence Distribution','NumberTitle','off');

subplot(1,2,1)
histogram(N_conv_all(valid), 'BinMethod','integers', 'FaceColor',[0.839 0.153 0.157]);
xlabel('Converged N_{segments}'); ylabel('Count');
title(sprintf('Distribution across %d feasible design points', sum(valid)));
grid on

subplot(1,2,2)
boxplot([N_conv_drag(valid), N_conv_length(valid), N_conv_dp(valid)], ...
        'Labels', {'drag\_tot','L\_solution','dp\_hx'});
ylabel('Converged N_{segments}');
title('Per-metric convergence N (which metric is slowest?)');
grid on

save('Nsegments_convergence_study.mat', 'design_points', 'N_conv_all', ...
     'N_conv_drag', 'N_conv_length', 'N_conv_dp', 'N_segments_candidates', 'tol', 'sample_status');

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
    fprintf('\n--- Part B: N_cool_seg convergence study at N_segments = %d (%d/%d) ---\n', ...
            N_seg_fixed, j, n_outer);

    all_entries_j = cell(n_samples,1);

    parfor i = 1:n_samples
        d2 = design_points.fr_dia(i);
        e  = design_points.air_side(i);
        r  = design_points.cool_side(i);

        collector = ResultsCollector();
        model_fun = @(Nc) capture( ...
            run_disc_model_fwdpass(use_DNS, N_seg_fixed, Nc, e, r, hx_theta, fan, fpr_init, ...
                M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
                d2, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, ...
                solve_T, false, false, tol_T), ...
            collector, e, r, d2, N_seg_fixed, Nc);

        try
            [N_conv, N_conv_per_metric, ~, status] = find_converged_Nsegments(model_fun, N_cool_seg_candidates, tol, metric_names);
            N_conv_coolseg_all(i,j)    = N_conv;
            N_conv_coolseg_drag(i,j)   = N_conv_per_metric.drag_tot;
            N_conv_coolseg_length(i,j) = N_conv_per_metric.L_solution;
            N_conv_coolseg_dp(i,j)     = N_conv_per_metric.dp_hx;
            sample_status_coolseg(i,j) = status;
        catch ME
            warning('Sample %d (N_segments=%d) failed entirely: %s', i, N_seg_fixed, ME.message);
            sample_status_coolseg(i,j) = "sample_error";
        end

        all_entries_j{i} = collector.entries;
    end

    % Save this j's model evaluations before moving to the next N_segments
    % value, same serial-after-parfor pattern as Part A.
    fprintf('Saving Part B (N_segments=%d) model evaluations to disk...\n', N_seg_fixed);
    n_saved_j = 0;
    for i = 1:n_samples
        for k = 1:numel(all_entries_j{i})
            ent = all_entries_j{i}{k};
            save_disc_results(ent.results, ent.e, ent.r, hx_theta, fan, fpr_init, ...
                M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
                ent.d2, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, ...
                solve_T, tol_T, ent.N_segments, ent.N_cool_seg, use_DNS, false);
            n_saved_j = n_saved_j + 1;
        end
    end
    fprintf('Saved %d model evaluations under results/ for N_segments=%d.\n', n_saved_j, N_seg_fixed);

    % Per-j failure reporting + feasibility diagnostic + distribution plot,
    % same structure as Part A.
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

    figure('Name', sprintf('Part B: Feasibility (N_segments=%d)', N_seg_fixed), 'NumberTitle', 'off');
    gscatter(design_points.fr_dia, design_points.cool_side, status_j);
    xlabel('fr\_dia'); ylabel('cool\_side');
    title(sprintf('Part B (N_{cool,seg} sweep, N_{segments}=%d): sample outcome', N_seg_fixed));

    if any(valid_j)
        fprintf('  Converged N_cool_seg (feasible subset, n=%d): median = %.1f, range %.0f-%.0f\n', ...
                sum(valid_j), median(N_conv_coolseg_all(valid_j,j)), ...
                min(N_conv_coolseg_all(valid_j,j)), max(N_conv_coolseg_all(valid_j,j)));

        figure('Name', sprintf('Part B: N_cool_seg Convergence Distribution (N_segments=%d)', N_seg_fixed), 'NumberTitle', 'off');

        subplot(1,2,1)
        histogram(N_conv_coolseg_all(valid_j,j), 'BinMethod','integers', 'FaceColor',[0.122 0.467 0.706]);
        xlabel('Converged N_{cool,seg}'); ylabel('Count');
        title(sprintf('N_{segments}=%d: distribution across %d feasible points', N_seg_fixed, sum(valid_j)));
        grid on

        subplot(1,2,2)
        boxplot([N_conv_coolseg_drag(valid_j,j), N_conv_coolseg_length(valid_j,j), N_conv_coolseg_dp(valid_j,j)], ...
                'Labels', {'drag\_tot','L\_solution','dp\_hx'});
        ylabel('Converged N_{cool,seg}');
        title('Per-metric convergence N (which metric is slowest?)');
        grid on
    else
        warning('Part B (N_segments=%d): no feasible/converged design points -- skipping distribution plot.', N_seg_fixed);
    end
end

%% ---- Cross-comparison: does converged N_cool_seg depend on N_segments? ---- %%
any_valid_coolseg = any(sample_status_coolseg == "ok", 1);
if any(any_valid_coolseg)
    figure('Name','Part B: N_cool_seg convergence vs. air-side resolution','NumberTitle','off');
    group_data = cell(1, n_outer);
    for j = 1:n_outer
        mask = sample_status_coolseg(:,j) == "ok";
        group_data{j} = N_conv_coolseg_all(mask, j);
    end
    % boxplot needs equal-length columns or the grouped (x,g) form -- use
    % the latter since each N_segments value can have a different number
    % of feasible/converged design points.
    x = vertcat(group_data{:});
    g = [];
    for j = 1:n_outer
        g = [g; repmat(N_segments_for_coolseg_study(j), numel(group_data{j}), 1)]; %#ok<AGROW>
    end
    boxplot(x, g);
    xlabel('N_{segments} (air-side resolution)'); ylabel('Converged N_{cool,seg}');
    title('Does the coolant-side resolution needed for convergence depend on the air-side resolution?');
    grid on
end

save('Ncoolseg_convergence_study.mat', 'design_points', 'N_conv_coolseg_all', ...
     'N_conv_coolseg_drag', 'N_conv_coolseg_length', 'N_conv_coolseg_dp', ...
     'N_cool_seg_candidates', 'N_segments_for_coolseg_study', 'tol', 'sample_status_coolseg');

%% ---- Local helper: tag a model evaluation into this iteration's collector ---- %%
% Passed through untouched as model_fun's return value -- find_converged_Nsegments
% only ever sees `res`, exactly as if this wrapper weren't here; the side
% effect (collector.add) is what makes every evaluation recoverable for
% the serial save pass afterward, see ResultsCollector.m.
function res = capture(res, collector, e, r, d2, N_segments, N_cool_seg)
collector.add(res, e, r, d2, N_segments, N_cool_seg);
end