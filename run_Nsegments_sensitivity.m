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

%% ---- N_segments candidates + convergence tolerance ---- %%
N_candidates = [1 2 5 10 15 20 30 40 50 75 100 150 200 250];
tol          = 0.01;  % 1%
metric_names = {'drag_tot', 'L_solution', 'dp_hx'};

%% ---- Monte Carlo sampling ---- %%
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

N_conv_all      = nan(n_samples,1);
N_conv_drag     = nan(n_samples,1);
N_conv_length   = nan(n_samples,1);
N_conv_dp       = nan(n_samples,1);
sample_status   = strings(n_samples,1); % "ok" | "no_convergence" | "infeasible_at_N" | "sample_error"

parfor i = 1:n_samples
    d2 = design_points.fr_dia(i);
    e  = design_points.air_side(i);
    r  = design_points.cool_side(i);

    model_fun = @(N) run_disc_model_fwdpass(false, N, e, r, hx_theta, fan, fpr_init, ...
        M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
        d2, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, tol_T);

    try
        [N_conv, N_conv_per_metric, ~, status] = find_converged_Nsegments(model_fun, N_candidates, tol, metric_names);
        N_conv_all(i)    = N_conv;
        N_conv_drag(i)   = N_conv_per_metric.drag_tot;
        N_conv_length(i) = N_conv_per_metric.L_solution;
        N_conv_dp(i)     = N_conv_per_metric.dp_hx;
        sample_status(i) = status;
    catch ME
        warning('Sample %d failed entirely (outside find_converged_Nsegments): %s', i, ME.message);
        sample_status(i) = "sample_error";
    end
end

%% ---- Failure reporting (before filtering to 'valid') ---- %%
n_ok             = sum(sample_status == "ok");
n_no_convergence = sum(sample_status == "no_convergence");
n_infeasible     = sum(sample_status == "infeasible_at_N");
n_sample_error   = sum(sample_status == "sample_error");

fprintf('\n=== Sample outcome breakdown (%d total design points) ===\n', n_samples);
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
figure('Name','Feasibility across design space','NumberTitle','off');
gscatter(design_points.fr_dia, design_points.air_side, sample_status);
xlabel('fr\_dia'); ylabel('air\_side');
title('Sample outcome across sampled design space');

%% ---- Summary statistics (feasible subset only) ---- %%
fprintf('\n=== N_segments convergence across %d FEASIBLE design points ===\n', sum(valid));
fprintf('Overall (all metrics simultaneous, tol=%.1f%%):\n', tol*100);
fprintf('  Min    = %.1f\n', min(N_conv_all(valid)));
fprintf('  Mean   = %.1f\n', mean(N_conv_all(valid)));
fprintf('  Median = %.1f\n', median(N_conv_all(valid)));
fprintf('  Max    = %.1f\n', max(N_conv_all(valid)));
fprintf('  Std    = %.1f\n', std(N_conv_all(valid)));

fprintf('\nPer-metric medians (diagnostic — which metric is slowest to converge):\n');
fprintf('  drag_tot    : median N = %.1f  (range %.0f-%.0f)\n', ...
        median(N_conv_drag(valid),'omitnan'), min(N_conv_drag(valid)), max(N_conv_drag(valid)));
fprintf('  L_solution  : median N = %.1f  (range %.0f-%.0f)\n', ...
        median(N_conv_length(valid),'omitnan'), min(N_conv_length(valid)), max(N_conv_length(valid)));
fprintf('  dp_hx       : median N = %.1f  (range %.0f-%.0f)\n', ...
        median(N_conv_dp(valid),'omitnan'), min(N_conv_dp(valid)), max(N_conv_dp(valid)));

%% ---- Plots ---- %%
figure('Name','N_segments Convergence Distribution','NumberTitle','off');

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
     'N_conv_drag', 'N_conv_length', 'N_conv_dp', 'N_candidates', 'tol', 'sample_status');