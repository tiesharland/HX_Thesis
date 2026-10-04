function plot_thermoturb_coverage(varargin)
%PLOT_THERMOTURB_COVERAGE  Visualize how a model run's DNS lookups sit
%relative to the thermoturb table's coverage.
%
%   plot_thermoturb_coverage()            -- 3 panels of 2D projections
%   plot_thermoturb_coverage('mode','3d') -- 1 panel, 3D scatter
%
%   Run this after a model run (anything that calls thermoturb_cached
%   with log_queries = true). It reads the accumulated log from
%   thermoturb_query_log('get') and the table's own laminar/turbulent
%   bounding boxes straight from thermoturb_table.mat, then plots every
%   logged query, colored by whether it fell inside the table's coverage
%   or was extrapolated/clamped in any of Re, Pr, Tb, Tw.
%
%   Call thermoturb_query_log('reset') before the run if you only want
%   that run's queries shown (otherwise the log accumulates across every
%   thermoturb_cached call made this session).

p = inputParser;
addParameter(p, 'mode', '2d', @(s) any(strcmpi(s, {'2d','3d'})));
parse(p, varargin{:});
mode = lower(p.Results.mode);

%% ---- Load table bounds directly ---- %%

table_file = fullfile(pwd, 'thermoturb_table.mat');
if ~isfile(table_file)
    error('plot_thermoturb_coverage: %s not found.', table_file);
end
loaded = load(table_file, 'tt');
tt     = loaded.tt;

sub_lam = tt(tt.regime == 0, :);
sub_tur = tt(tt.regime == 1, :);

box_lam = struct('Re', [min(sub_lam.Re), max(sub_lam.Re)], ...
                  'Pr', [min(sub_lam.Pr), max(sub_lam.Pr)], ...
                  'Tb', [min(sub_lam.Tb), max(sub_lam.Tb)], ...
                  'Tw', [min(sub_lam.Tw), max(sub_lam.Tw)]);
box_tur = struct('Re', [min(sub_tur.Re), max(sub_tur.Re)], ...
                  'Pr', [min(sub_tur.Pr), max(sub_tur.Pr)], ...
                  'Tb', [min(sub_tur.Tb), max(sub_tur.Tb)], ...
                  'Tw', [min(sub_tur.Tw), max(sub_tur.Tw)]);

%% ---- Load logged queries ---- %%

log_tbl = thermoturb_query_log('get');
if isempty(log_tbl) || height(log_tbl) == 0
    error(['plot_thermoturb_coverage: the query log is empty. Run the model ' ...
           'first with thermoturb_cached''s log_queries flag set to true.']);
end

Re = log_tbl.Re; Tb = log_tbl.Tb; Tw = log_tbl.Tw;
any_extrap = log_tbl.any_extrap;

fprintf('plot_thermoturb_coverage: %d logged queries, %d (%.1f%%) extrapolated/clamped.\n', ...
        height(log_tbl), sum(any_extrap), 100*mean(any_extrap));

%% ---- Colours ---- %%

col_lam    = [0.122 0.471 0.706];   % blue   — laminar table box
col_tur    = [0.173 0.627 0.173];   % green  — turbulent table box
col_in     = [0.3 0.3 0.3];         % gray   — in-coverage query
col_out    = [0.839 0.153 0.157];   % red    — extrapolated/clamped query

in_mask  = ~any_extrap;
out_mask =  any_extrap;

%% ---- Plot ---- %%

if strcmp(mode, '2d')

    figure('Name', 'thermoturb DNS coverage', 'NumberTitle', 'off');

    panels = {
        'Re', 'Tb', 'Re_{bulk}', 'T_b [K]';
        'Re', 'Tw', 'Re_{bulk}', 'T_w [K]';
        'Tb', 'Tw', 'T_b [K]',   'T_w [K]'
    };

    for k = 1:3
        subplot(1, 3, k)
        hold on

        xf = panels{k,1}; yf = panels{k,2};
        x_all = log_tbl.(xf); y_all = log_tbl.(yf);

        draw_box(box_lam.(xf), box_lam.(yf), col_lam, 'Laminar table');
        draw_box(box_tur.(xf), box_tur.(yf), col_tur, 'Turbulent table');

        scatter(x_all(in_mask),  y_all(in_mask),  18, col_in,  'filled', ...
                'DisplayName', 'In coverage', 'MarkerFaceAlpha', 0.5)
        scatter(x_all(out_mask), y_all(out_mask), 18, col_out, 'filled', ...
                'DisplayName', 'Extrapolated/clamped')

        xlabel(panels{k,3}); ylabel(panels{k,4});
        title(sprintf('%s vs %s', panels{k,4}, panels{k,3}))
        grid on
        if k == 1, legend('Location', 'best', 'FontSize', 7); end
    end

    sgtitle('DNS query coverage vs. thermoturb table bounds')

else % 3d

    figure('Name', 'thermoturb DNS coverage (3D)', 'NumberTitle', 'off');
    hold on

    draw_cuboid(box_lam.Re, box_lam.Tb, box_lam.Tw, col_lam, 'Laminar table');
    draw_cuboid(box_tur.Re, box_tur.Tb, box_tur.Tw, col_tur, 'Turbulent table');

    scatter3(Re(in_mask),  Tb(in_mask),  Tw(in_mask),  18, col_in,  'filled', ...
             'DisplayName', 'In coverage', 'MarkerFaceAlpha', 0.5)
    scatter3(Re(out_mask), Tb(out_mask), Tw(out_mask), 18, col_out, 'filled', ...
             'DisplayName', 'Extrapolated/clamped')

    xlabel('Re_{bulk}'); ylabel('T_b [K]'); zlabel('T_w [K]');
    title('DNS query coverage vs. thermoturb table bounds')
    legend('Location', 'best')
    grid on
    view(3)

end

end

%% ---- Local helpers ---- %%

function draw_box(xrng, yrng, col, name)
rectangle('Position', [xrng(1), yrng(1), diff(xrng), diff(yrng)], ...
          'EdgeColor', col, 'LineWidth', 1.5, 'LineStyle', '--');
% Invisible proxy line so the box gets a legend entry.
plot(nan, nan, '--', 'Color', col, 'LineWidth', 1.5, 'DisplayName', name);
end

function draw_cuboid(xrng, yrng, zrng, col, name)
[X, Y, Z] = ndgrid(xrng, yrng, zrng);
verts = [X(:), Y(:), Z(:)];
edges = [1 2; 3 4; 5 6; 7 8; ...   % along x
         1 3; 2 4; 5 7; 6 8; ...   % along y
         1 5; 2 6; 3 7; 4 8];      % along z
for i = 1:size(edges, 1)
    v = verts(edges(i,:), :);
    plot3(v(:,1), v(:,2), v(:,3), '--', 'Color', col, 'LineWidth', 1.2, ...
          'HandleVisibility', 'off');
end
% Invisible proxy for the legend.
plot3(nan, nan, nan, '--', 'Color', col, 'LineWidth', 1.5, 'DisplayName', name);
end