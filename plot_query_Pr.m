function plot_query_Pr(log_tbl)
%PLOT_QUERY_PR  Prandtl number of every logged thermoturb query.
%
%   plot_query_Pr()          -- live log from thermoturb_query_log('get')
%   plot_query_Pr(log_tbl)   -- a saved log, e.g. results.thermoturb_log
%
% Left: Pr per query in chronological order. Right: histogram. Dashed lines
% mark the Pr values present in thermoturb_table.mat (if found in the
% current folder), so you can see whether the queries sit between them.

if nargin < 1 || isempty(log_tbl)
    log_tbl = thermoturb_query_log('get');
end
if isempty(log_tbl) || height(log_tbl) == 0
    error('plot_query_Pr: the query log is empty.');
end

Pr = log_tbl.Pr;

fprintf('%d queries: Pr min = %.4f, median = %.4f, max = %.4f\n', ...
    numel(Pr), min(Pr), median(Pr), max(Pr));

Pr_grid = [];
table_file = fullfile(pwd, 'thermoturb_table.mat');
if isfile(table_file)
    loaded  = load(table_file, 'tt');
    Pr_grid = unique(round(loaded.tt.Pr, 4));
end

figure('Name', 'Prandtl number of logged queries', 'NumberTitle', 'off');

subplot(1, 2, 1)
plot(Pr, '.', 'MarkerSize', 6); hold on
for g = Pr_grid'
    yline(g, ':', 'Color', [0.6 0.6 0.6]);
end
xlabel('Query #'); ylabel('Pr'); grid on
title('Pr per query')

subplot(1, 2, 2)
histogram(Pr, 30); hold on
for g = Pr_grid'
    xline(g, ':', 'Color', [0.6 0.6 0.6]);
end
xlabel('Pr'); ylabel('Count'); grid on
title('Distribution')

end