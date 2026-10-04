function plot_thermoturb_coverage(varargin)
%PLOT_THERMOTURB_COVERAGE  Visualize how a model run's DNS lookups sit
%relative to the thermoturb table's coverage.
%
%   plot_thermoturb_coverage()                 -- 3 panels of 2D projections
%   plot_thermoturb_coverage('mode','3d')      -- 1 panel, 3D scatter + wireframe box
%   plot_thermoturb_coverage('highlight', K)   -- color the LAST K logged
%                                                  queries bright red
%                                                  (e.g. K = N_segments*N_seg_cool
%                                                  for the final design
%                                                  iteration of a sweep),
%                                                  drawn on top of the rest.
%                                                  K = 0 (default): off,
%                                                  no query stands out.
%                                                  Either way the bulk
%                                                  points are colored by
%                                                  chronological order on
%                                                  MATLAB's default colormap.
%   plot_thermoturb_coverage('log', log_tbl)   -- plot from a previously
%                                                  saved query log (e.g.
%                                                  results.thermoturb_log
%                                                  from a saved results
%                                                  file) instead of the
%                                                  live persistent log.
%                                                  Lets you replot any old
%                                                  run's coverage without
%                                                  re-running the model.
%
%   Run this after a model run (anything that calls thermoturb_cached
%   with log_queries = true). By default it reads the accumulated log
%   from thermoturb_query_log('get'); pass 'log' to use a saved one
%   instead. Either way, the table's own coverage bounds are read
%   straight from thermoturb_table.mat.
%
%   "Coverage" here means: a query's (Pr,Tb,Tw) cell is inside the
%   sampled grid. Re is NOT part of that test -- within an in-grid
%   (Pr,Tb,Tw) combo, the whole laminar -> transition-blend ->
%   turbulent-extrapolated -> turbulent-DNS path is considered covered,
%   since thermoturb_cached always returns a value there (it's only
%   extrapolating in Re, same combo). Only a (Pr,Tb,Tw) that falls
%   outside the sampled grid entirely is "not covered".
%
%   Coverage box: not every (Pr,Tb,Tw) combo's Re sweep spans the same
%   range, and Pr varies the range only slightly, so rather than drawing
%   one flat min/max rectangle, the box is built from the 4 corners of
%   the Tb-Tw grid (min/max Tb x min/max Tw). At each corner, across the
%   Pr values sampled there, the Pr-slice with the NARROWEST Re span is
%   used as that corner's representative range (the conservative choice
%   given Pr's small effect). That gives three Re values per corner:
%   laminar start, where real turbulent DNS data ends, and turbulent end
%   (incl. extrapolated). Re=2300/3000 (the fixed laminar/turbulent
%   demarcation) are drawn as plain reference lines.
%
%   In 2D, each boundary runs as a straight line between corner values;
%   a 2D panel can only show Re against one of Tb/Tw at a time, so the
%   OTHER is collapsed by drawing the boundary at its two grid extremes
%   as two separate lines rather than merging them.
%
%   In 3D, the three boundaries form one wireframe box: vertical edges
%   run laminar-start -> DNS-end -> turbulent-end at each of the 4
%   corners, with a rectangular ring at each of those three levels. Only
%   the top (turbulent-end) ring is also filled in as a light gray plane
%   -- the bottom and middle levels stay wireframe-only.
%
%   Call thermoturb_query_log('reset') before the run if you only want
%   that run's queries shown (otherwise the log accumulates across every
%   thermoturb_cached call made this session -- including, e.g., a
%   separate Re-sweep test script run earlier).

p = inputParser;
addParameter(p, 'mode', '2d', @(s) any(strcmpi(s, {'2d','3d'})));
addParameter(p, 'highlight', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'log', [], @(t) isempty(t) || istable(t));
parse(p, varargin{:});
mode      = lower(p.Results.mode);
highlight = round(p.Results.highlight);
log_input = p.Results.log;

% Fixed laminar/turbulent demarcation used by thermoturb_cached -- only
% used here to draw reference lines, not recomputed.
Re_lam_max = 2300;
Re_tur_min = 3000;

%% ---- Load table and compute the 4-corner coverage boundaries ---- %%

table_file = fullfile(pwd, 'thermoturb_table.mat');
if ~isfile(table_file)
    error('plot_thermoturb_coverage: %s not found.', table_file);
end
loaded = load(table_file, 'tt');
tt     = loaded.tt;

grid_TbTw = unique(tt{:, {'Tb','Tw'}}, 'rows');
Tb_min = min(grid_TbTw(:,1)); Tb_max = max(grid_TbTw(:,1));
Tw_min = min(grid_TbTw(:,2)); Tw_max = max(grid_TbTw(:,2));

% Corner order: 1=(Tb_min,Tw_min) 2=(Tb_min,Tw_max) 3=(Tb_max,Tw_min) 4=(Tb_max,Tw_max)
corners = [Tb_min Tw_min; Tb_min Tw_max; Tb_max Tw_min; Tb_max Tw_max];

lam_lo = zeros(4,1);   % laminar start, per corner (box's lower edge)
tur_hi = zeros(4,1);   % turbulent end (incl. extrapolated), per corner (box's upper edge)
tur_lo = zeros(4,1);   % where real turbulent DNS data ends, per corner (middle level)
for i = 1:4
    [lo, ~]    = corner_re_range(tt, corners(i,1), corners(i,2), 0);
    [lo2, hi2] = corner_re_range(tt, corners(i,1), corners(i,2), 1);
    lam_lo(i) = lo;
    tur_lo(i) = lo2;
    tur_hi(i) = hi2;
end

%% ---- Load logged queries ---- %%

if isempty(log_input)
    log_tbl = thermoturb_query_log('get');
else
    log_tbl = log_input;
end
if isempty(log_tbl) || height(log_tbl) == 0
    error(['plot_thermoturb_coverage: the query log is empty. Run the model ' ...
           'first with thermoturb_cached''s log_queries flag set to true, or ' ...
           'pass a saved log via the ''log'' argument (e.g. results.thermoturb_log).']);
end

n_pts = height(log_tbl);

% "Not covered" = the query's (Pr,Tb,Tw) cell fell outside the sampled
% grid. Re_extrap does NOT count -- extrapolating in Re within an
% in-grid combo is covered.
not_covered = log_tbl.Pr_extrap | log_tbl.Tb_extrap | log_tbl.Tw_extrap;

fprintf('plot_thermoturb_coverage: %d logged queries, %d (%.1f%%) outside the (Pr,Tb,Tw) grid.\n', ...
        n_pts, sum(not_covered), 100*mean(not_covered));

% Chronological order color, with the trailing 'highlight' queries (e.g.
% the final design iteration's N_segments*N_seg_cool calls) singled out
% in bright red and drawn on top. The bulk points always stay on
% MATLAB's default colormap, whether or not highlighting is on.
highlight    = min(highlight, n_pts);
is_highlight = false(n_pts, 1);
if highlight > 0
    is_highlight(end-highlight+1:end) = true;
end
bulk_mask     = ~is_highlight;
order_idx     = (1:n_pts)';
bulk_cmap     = parula(256);   % MATLAB default, always
col_highlight = [1 0 0];

%% ---- Colours ---- %%

col_box    = [0.55 0.55 0.55];   % box edges (laminar start / turbulent end)
col_dnsend = [0.55 0.55 0.55];   % turbulent-DNS-ends reference line
col_thresh = [0.55 0.55 0.55];   % Re=2300 / Re=3000 reference lines
col_grid   = [0.55 0.55 0.55];   % Tb-Tw sampled-grid footprint

%% ---- Plot ---- %%

if strcmp(mode, '2d')

    fig = figure('Name', 'thermoturb DNS coverage', 'NumberTitle', 'off');

    panels = {
        'Re', 'Tb', 'Re_{bulk}', 'T_b [K]';
        'Re', 'Tw', 'Re_{bulk}', 'T_w [K]';
        'Tb', 'Tw', 'T_b [K]',   'T_w [K]'
    };

    for k = 1:3
        ax = subplot(1, 3, k);
        hold(ax, 'on')

        xf = panels{k,1}; yf = panels{k,2};
        x_all = log_tbl.(xf); y_all = log_tbl.(yf);

        % Bulk points first (so reference lines drawn after stay visible
        % on top of dense clusters), colormap set per-axes right after.
        scatter(ax, x_all(bulk_mask), y_all(bulk_mask), 12, order_idx(bulk_mask), 'o', 'filled', ...
                'MarkerFaceAlpha', 0.35, 'HandleVisibility', 'off');
        colormap(ax, bulk_cmap);

        if strcmp(xf, 'Re') && strcmp(yf, 'Tb')
            draw_corner_lines(ax, Tb_min, Tb_max, lam_lo([1 3]), lam_lo([2 4]), col_box,    '-',  1.5, 'Coverage start (laminar)');
            draw_corner_lines(ax, Tb_min, Tb_max, tur_hi([1 3]), tur_hi([2 4]), col_box,    '-',  1.5, 'Coverage end (turbulent, incl. extrap.)');
            draw_corner_lines(ax, Tb_min, Tb_max, tur_lo([1 3]), tur_lo([2 4]), col_dnsend, '--', 1.0, 'Turbulent DNS data ends');
            xline(ax, Re_lam_max, ':', 'Color', col_thresh, 'HandleVisibility', 'off');
            xline(ax, Re_tur_min, ':', 'Color', col_thresh, 'HandleVisibility', 'off');
        elseif strcmp(xf, 'Re') && strcmp(yf, 'Tw')
            draw_corner_lines(ax, Tw_min, Tw_max, lam_lo([1 2]), lam_lo([3 4]), col_box,    '-',  1.5, 'Coverage start (laminar)');
            draw_corner_lines(ax, Tw_min, Tw_max, tur_hi([1 2]), tur_hi([3 4]), col_box,    '-',  1.5, 'Coverage end (turbulent, incl. extrap.)');
            draw_corner_lines(ax, Tw_min, Tw_max, tur_lo([1 2]), tur_lo([3 4]), col_dnsend, '--', 1.0, 'Turbulent DNS data ends');
            xline(ax, Re_lam_max, ':', 'Color', col_thresh, 'HandleVisibility', 'off');
            xline(ax, Re_tur_min, ':', 'Color', col_thresh, 'HandleVisibility', 'off');
        else
            rectangle(ax, 'Position', [Tb_min, Tw_min, Tb_max-Tb_min, Tw_max-Tw_min], ...
                      'EdgeColor', col_grid, 'LineStyle', '--', 'LineWidth', 1.2);
            plot(ax, nan, nan, '--', 'Color', col_grid, 'LineWidth', 1.2, ...
                 'DisplayName', 'Sampled (Tb,Tw) grid');
        end

        % Highlighted points drawn last, on top of everything.
        if any(is_highlight)
            scatter(ax, x_all(is_highlight), y_all(is_highlight), 30, col_highlight, 'o', 'filled', ...
                    'MarkerEdgeColor', 'k', 'DisplayName', 'Final design iteration', ...
                    'HandleVisibility', tern(k==1));
        end

        xlabel(ax, panels{k,3}); ylabel(ax, panels{k,4});
        title(ax, sprintf('%s vs %s', panels{k,4}, panels{k,3}))
        grid(ax, 'on')
        cb = colorbar(ax);
        cb.Label.String = 'Query order (chronological)';
        if k == 1, legend(ax, 'Location', 'best', 'FontSize', 7); end
    end

    sgtitle(fig, 'DNS query coverage vs. thermoturb table bounds')

else % 3d

    fig = figure('Name', 'thermoturb DNS coverage (3D)', 'NumberTitle', 'off');
    ax = axes(fig);
    hold(ax, 'on')

    scatter3(ax, log_tbl.Re(bulk_mask), log_tbl.Tb(bulk_mask), log_tbl.Tw(bulk_mask), ...
             12, order_idx(bulk_mask), 'o', 'filled', 'MarkerFaceAlpha', 0.35, 'HandleVisibility', 'off');
    colormap(ax, bulk_cmap);

    draw_wireframe_box(ax, lam_lo, tur_lo, tur_hi, corners, col_box, 'Coverage box (laminar start / DNS end / turbulent end)');

    if any(is_highlight)
        scatter3(ax, log_tbl.Re(is_highlight), log_tbl.Tb(is_highlight), log_tbl.Tw(is_highlight), ...
                 36, col_highlight, 'o', 'filled', 'MarkerEdgeColor', 'k', ...
                 'DisplayName', 'Final design iteration');
    end

    xlabel(ax, 'Re_{bulk}'); ylabel(ax, 'T_b [K]'); zlabel(ax, 'T_w [K]');
    title(ax, 'DNS query coverage vs. thermoturb table bounds')
    cb = colorbar(ax);
    cb.Label.String = 'Query order (chronological)';
    legend(ax, 'Location', 'best')
    grid(ax, 'on')
    view(ax, 3)

end

end

%% ---- Local helper: narrowest-Pr-span Re range at one (Tb,Tw) corner ---- %%
function [re_lo, re_hi] = corner_re_range(tt, Tb_val, Tw_val, regime)
mask = tt.regime == regime & abs(tt.Tb - Tb_val) < 0.1 & abs(tt.Tw - Tw_val) < 0.1;
sub  = tt(mask, :);
Pr_vals = unique(sub.Pr);
n = numel(Pr_vals);
los = zeros(n,1); his = zeros(n,1); spans = zeros(n,1);
for j = 1:n
    rows   = abs(sub.Pr - Pr_vals(j)) < 1e-4;
    los(j) = min(sub.Re(rows));
    his(j) = max(sub.Re(rows));
    spans(j) = his(j) - los(j);
end
[~, sel] = min(spans);
re_lo = los(sel);
re_hi = his(sel);
end

%% ---- Local helper: two straight boundary lines (one per hidden-axis extreme) ---- %%
% y_lo/y_hi are the two grid extremes of the axis NOT shown as Re (Tb or
% Tw); vals_fixed_A/B are the two corner Re-values spanning y_lo to y_hi
% at each of the other (hidden) axis's two extremes.
function draw_corner_lines(ax, y_lo, y_hi, vals_fixed_A, vals_fixed_B, col, ln, lw, name)
plot(ax, vals_fixed_A, [y_lo y_hi], ln, 'Color', col, 'LineWidth', lw, 'HandleVisibility', 'off');
plot(ax, vals_fixed_B, [y_lo y_hi], ln, 'Color', col, 'LineWidth', lw, 'HandleVisibility', 'off');
% Single legend proxy for the pair.
plot(ax, nan, nan, ln, 'Color', col, 'LineWidth', lw, 'DisplayName', name);
end

%% ---- Local helper: 3D wireframe box through the 3 corner levels (3D) ---- %%
% corners rows: 1=(Tb_min,Tw_min) 2=(Tb_min,Tw_max) 3=(Tb_max,Tw_min) 4=(Tb_max,Tw_max)
% lo/mid/hi are each 4x1, the Re value of that level at each corner.
% Vertical edges run lo -> mid -> hi at each corner; each level also gets
% a rectangular ring around its 4 corners. Only the top (hi) ring is also
% filled in as a translucent plane.
function draw_wireframe_box(ax, lo, mid, hi, corners, col, name)
ring_order = [1 3 4 2 1];   % closed loop around the rectangle perimeter

plot3(ax, lo(ring_order),  corners(ring_order,1), corners(ring_order,2), '-',  'Color', col, 'LineWidth', 1.2, 'HandleVisibility', 'off');
plot3(ax, mid(ring_order), corners(ring_order,1), corners(ring_order,2), '--', 'Color', col, 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot3(ax, hi(ring_order),  corners(ring_order,1), corners(ring_order,2), '-',  'Color', col, 'LineWidth', 1.2, 'HandleVisibility', 'off');

for i = 1:4
    plot3(ax, [lo(i) mid(i) hi(i)], corners(i,1)*[1 1 1], corners(i,2)*[1 1 1], '-', ...
          'Color', col, 'LineWidth', 1.0, 'HandleVisibility', 'off');
end

patch(ax, 'XData', hi(ring_order(1:4)), 'YData', corners(ring_order(1:4),1), 'ZData', corners(ring_order(1:4),2), ...
      'FaceColor', col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

plot3(ax, nan, nan, nan, '-', 'Color', col, 'LineWidth', 1.5, 'DisplayName', name);
end

%% ---- Local helper: 'on'/'off' for HandleVisibility from a logical ---- %%
function s = tern(tf)
if tf, s = 'on'; else, s = 'off'; end
end