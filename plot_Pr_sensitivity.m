function plot_Pr_sensitivity(log_tbl, n_sample)
%PLOT_PR_SENSITIVITY  How much does the heat transfer coefficient change
%with Prandtl number, at the conditions the model actually queried?
%
%   plot_Pr_sensitivity()               -- live log, 300 sampled queries
%   plot_Pr_sensitivity(log_tbl)        -- saved log, e.g. results.thermoturb_log
%   plot_Pr_sensitivity(log_tbl, 1000)  -- sample more queries
%
% For each sampled logged query (Re, Tb, Tw held fixed) this looks up Nu from
% thermoturb_table.mat across the table's whole Pr range. Since h = Nu*k/d_h
% and neither k nor d_h depends on the table's Pr, the % change in Nu is the
% % change in h.
%
% Left:  % change in h vs Pr, relative to Nu at the median logged Pr, one
%        line per sampled query. Dashed lines mark the min/max Pr actually
%        logged; dots mark each query's own Pr.
% Right: histogram, over the sampled queries, of the % change in h between
%        the smallest and largest logged Pr.
%
% This reimplements the lookup of thermoturb_cached (Re-interpolation within
% each (Pr,Tb,Tw) group, then linear interpolation over Pr/Tb/Tw, inputs
% clamped to the table) instead of calling it, because every call to
% thermoturb_cached would append to the live query log. Run it on the current
% table, before extending it: a table with only Pr 0.70/0.71 in part of the
% (Tb,Tw) region would not show the real Pr dependence there.

if nargin < 1 || isempty(log_tbl)
    log_tbl = thermoturb_query_log('get');
end
if nargin < 2 || isempty(n_sample)
    n_sample = 300;
end
if isempty(log_tbl) || height(log_tbl) == 0
    error('plot_Pr_sensitivity: the query log is empty.');
end

table_file = fullfile(pwd, 'thermoturb_table.mat');
if ~isfile(table_file)
    error('plot_Pr_sensitivity: %s not found.', table_file);
end
loaded = load(table_file, 'tt');
tt     = loaded.tt;

Re_lam_max = 2300;
Re_tur_min = 3000;

Gl = build_groups(tt(tt.regime == 0, :));
Gt = build_groups(tt(tt.regime == 1, :));

Pr_tab = [min([Gl.Pr_rng(1) Gt.Pr_rng(1)]), max([Gl.Pr_rng(2) Gt.Pr_rng(2)])];
Pr_log = log_tbl.Pr;
Pr_lo  = max(min(Pr_log), Pr_tab(1));
Pr_hi  = min(max(Pr_log), Pr_tab(2));
Pr_ref = min(max(median(Pr_log), Pr_tab(1)), Pr_tab(2));

fprintf('Logged Pr: min %.4f, median %.4f, max %.4f (table covers %.3f-%.3f)\n', ...
        min(Pr_log), median(Pr_log), max(Pr_log), Pr_tab(1), Pr_tab(2));

Pr_eval = linspace(Pr_tab(1), Pr_tab(2), 41)';
Pr_vec  = [Pr_eval; Pr_ref; Pr_lo; Pr_hi];

rng(1);
n_pts    = height(log_tbl);
sel      = randperm(n_pts, min(n_sample, n_pts));
n_sel    = numel(sel);
rel_curve = nan(n_sel, numel(Pr_eval));
pct_range = nan(n_sel, 1);
pct_table = nan(n_sel, 1);
dot_y     = nan(n_sel, 1);

for k = 1:n_sel
    q  = sel(k);
    Nu = nu_query(Gl, Gt, log_tbl.Re(q), Pr_vec, log_tbl.Tb(q), log_tbl.Tw(q), Re_lam_max, Re_tur_min);

    Nu_curve = Nu(1:numel(Pr_eval));
    Nu_ref   = Nu(numel(Pr_eval) + 1);
    Nu_lo    = Nu(numel(Pr_eval) + 2);
    Nu_hi    = Nu(numel(Pr_eval) + 3);

    rel_curve(k, :) = 100 * (Nu_curve / Nu_ref - 1);
    pct_range(k)    = 100 * (Nu_hi - Nu_lo) / Nu_lo;
    pct_table(k)    = 100 * (Nu_curve(end) - Nu_curve(1)) / Nu_curve(1);

    Pr_q   = min(max(log_tbl.Pr(q), Pr_tab(1)), Pr_tab(2));
    dot_y(k) = interp1(Pr_eval, rel_curve(k, :), Pr_q);
end

fprintf('Change in h between logged Pr min and max (%d sampled queries):\n', n_sel);
fprintf('  median %.3f %%,  95th pct |.| %.3f %%,  max |.| %.3f %%\n', ...
        median(pct_range, 'omitnan'), prctile(abs(pct_range), 95), max(abs(pct_range)));
fprintf('Change in h across the whole table Pr range (%.3f -> %.3f):\n', Pr_tab(1), Pr_tab(2));
fprintf('  median %.3f %%,  95th pct |.| %.3f %%,  max |.| %.3f %%\n', ...
        median(pct_table, 'omitnan'), prctile(abs(pct_table), 95), max(abs(pct_table)));

figure('Name', 'h sensitivity to Prandtl number', 'NumberTitle', 'off');

subplot(1, 2, 1)
plot(Pr_eval, rel_curve', '-', 'Color', [0.75 0.75 0.75]); hold on
plot(Pr_log(sel), dot_y, '.', 'Color', [0.839 0.153 0.157], 'MarkerSize', 8);
xline(Pr_lo, '--k'); xline(Pr_hi, '--k');
xlabel('Pr'); ylabel('Change in h vs. median logged Pr [%]'); grid on
title('h vs. Pr at each query''s own Re, T_b, T_w')

subplot(1, 2, 2)
histogram(pct_range, 30, 'FaceColor', [0.839 0.153 0.157]);
xlabel(sprintf('Change in h, Pr %.3f -> %.3f [%%]', Pr_lo, Pr_hi)); ylabel('Count'); grid on
title('Spread over the Pr range actually queried')

end

%% ---- Local helper: per-(Pr,Tb,Tw) groups of one regime, Re-sorted ---- %%
function G = build_groups(sub)
[grps, ~, ic] = unique(sub{:, {'Pr','Tb','Tw'}}, 'rows');
[~, ord] = sortrows([ic, sub.Re]);
ic = ic(ord); Re = sub.Re(ord); Nu = sub.Nu(ord);
starts = [1; find(diff(ic) ~= 0) + 1];
stops  = [starts(2:end) - 1; numel(ic)];

n = size(grps, 1);
G.Pr = grps(:,1); G.Tb = grps(:,2); G.Tw = grps(:,3);
G.Re = cell(n, 1); G.Nu = cell(n, 1);
for i = 1:n
    r = Re(starts(i):stops(i));
    v = Nu(starts(i):stops(i));
    [r, iu] = unique(r);
    G.Re{i} = r;
    G.Nu{i} = v(iu);
end
G.Pr_rng = [min(G.Pr) max(G.Pr)];
G.Tb_rng = [min(G.Tb) max(G.Tb)];
G.Tw_rng = [min(G.Tw) max(G.Tw)];
G.F = scatteredInterpolant(G.Pr, G.Tb, G.Tw, zeros(n, 1), 'linear', 'nearest');
end

%% ---- Local helper: Nu for one regime at several Pr, fixed Re/Tb/Tw ---- %%
function Nu = nu_at(G, Re, Pr_vec, Tb, Tw)
n = numel(G.Pr);
Nu_g = zeros(n, 1);
for i = 1:n
    Nu_g(i) = interp1(G.Re{i}, G.Nu{i}, Re, 'linear', 'extrap');
end
G.F.Values = Nu_g;
Pr_q = min(max(Pr_vec, G.Pr_rng(1)), G.Pr_rng(2));
Tb_q = min(max(Tb, G.Tb_rng(1)), G.Tb_rng(2));
Tw_q = min(max(Tw, G.Tw_rng(1)), G.Tw_rng(2));
Nu = G.F(Pr_q, Tb_q * ones(size(Pr_q)), Tw_q * ones(size(Pr_q)));
end

%% ---- Local helper: same laminar / blend / turbulent split as thermoturb_cached ---- %%
function Nu = nu_query(Gl, Gt, Re, Pr_vec, Tb, Tw, Re_lam_max, Re_tur_min)
if Re <= Re_lam_max
    Nu = nu_at(Gl, Re, Pr_vec, Tb, Tw);
elseif Re >= Re_tur_min
    Nu = nu_at(Gt, Re, Pr_vec, Tb, Tw);
else
    w  = 0.5 * (1 - cos(pi * (Re - Re_lam_max) / (Re_tur_min - Re_lam_max)));
    Nu = (1 - w) * nu_at(Gl, Re, Pr_vec, Tb, Tw) + w * nu_at(Gt, Re, Pr_vec, Tb, Tw);
end
end