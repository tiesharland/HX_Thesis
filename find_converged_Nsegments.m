function [N_conv, N_conv_per_metric, results_table, status] = find_converged_Nsegments(model_fun, N_candidates, tol, metric_names)
% FIND_CONVERGED_NSEGMENTS  Find the smallest N_segments where relative
% change in tracked metrics vs. the previous N_segments falls below tol.
%
% If model_fun(N) fails at ANY candidate N, evaluation stops immediately
% (higher N is assumed equally infeasible, since N_segments is a numerical
% resolution parameter and does not change the underlying design's
% feasibility). The design point is flagged for discarding via `status`.
%
% Inputs:
%   model_fun     - function handle, model_fun(N) returns a results struct
%   N_candidates  - increasing vector of integer N_segments to test
%   tol           - relative change tolerance, e.g. 0.01 for 1%
%   metric_names  - cellstr of struct fieldnames to track
%
% Outputs:
%   N_conv             - smallest N where ALL metrics converge (NaN if not found)
%   N_conv_per_metric  - struct, per-metric convergence N (NaN if not found)
%   results_table       - table of tested N_segments, metric values, rel. change
%                         (only contains rows up to the point of failure, if any)
%   status              - "ok" | "no_convergence" | "infeasible_at_N"
%                         "infeasible_at_N" means model_fun failed at some N
%                         and evaluation was stopped early

n_metrics = numel(metric_names);
n_N = numel(N_candidates);

metric_vals = nan(n_N, n_metrics);
status = "ok";
N_failed_at = NaN;

for i = 1:n_N
    N = N_candidates(i);
    try
        res = model_fun(N);
        for m = 1:n_metrics
            metric_vals(i,m) = res.(metric_names{m});
        end
    catch ME
        warning('Model evaluation failed at N_segments=%d — assuming design is infeasible, skipping higher N. (%s)', ...
                 N, ME.message);
        status = "infeasible_at_N";
        N_failed_at = N;
        break  % do NOT try higher N for this design point
    end
end

% Trim to only the rows actually evaluated (relevant if we broke out early)
n_evaluated = find(~all(isnan(metric_vals),2), 1, 'last');
if isempty(n_evaluated), n_evaluated = 0; end
N_candidates_eval = N_candidates(1:n_evaluated);
metric_vals_eval  = metric_vals(1:n_evaluated,:);

if status == "infeasible_at_N"
    % Design point discarded entirely — no point computing convergence off a partial sweep
    N_conv = NaN;
    N_conv_per_metric = cell2struct(num2cell(nan(n_metrics,1)), metric_names, 1);
    results_table = table(N_candidates_eval(:), metric_vals_eval, nan(numel(N_candidates_eval), n_metrics), ...
        'VariableNames', {'N_segments', 'metric_values', 'relative_change'});
    results_table.Properties.UserData = struct('failed_at_N', N_failed_at);
    return
end

% --- Normal path: full sweep succeeded, compute convergence as before ---
rel_change = nan(n_N, n_metrics);
for m = 1:n_metrics
    for i = 2:n_N
        prev = metric_vals(i-1, m);
        curr = metric_vals(i, m);
        if abs(prev) > eps
            rel_change(i,m) = abs(curr - prev) / abs(prev);
        else
            rel_change(i,m) = abs(curr - prev);
        end
    end
end

N_conv_per_metric = struct();
for m = 1:n_metrics
    idx = find(rel_change(:,m) < tol, 1, 'first');
    if isempty(idx)
        N_conv_per_metric.(metric_names{m}) = NaN;
    else
        N_conv_per_metric.(metric_names{m}) = N_candidates(idx);
    end
end

all_converged = all(rel_change < tol, 2);
idx_all = find(all_converged, 1, 'first');
if isempty(idx_all)
    N_conv = NaN;
    status = "no_convergence";
    warning('No N_segments in the candidate list satisfied the %.1f%% tolerance for all metrics.', tol*100);
else
    N_conv = N_candidates(idx_all);
end

results_table = table(N_candidates(:), metric_vals, rel_change, ...
    'VariableNames', {'N_segments', 'metric_values', 'relative_change'});

end