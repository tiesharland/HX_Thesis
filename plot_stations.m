function figs_out = plot_stations(results, model_type, varargin)
% plot_stations  Add one model's TMS station data to a shared set of
% comparison figures - same calling convention as plot_HX.
%
% Usage (one call per dataset, same figs struct passed back in each time):
%   figs = plot_stations(results_D,    'discretised');
%   figs = plot_stations(results_L,    'lumped',      'figs', figs);
%   figs = plot_stations(results_DNS,  'dns',         'figs', figs);
%   figs = plot_stations(results_2D,   '2D',          'figs', figs);
%   figs = plot_stations(results_2DNS, '2DNS',        'figs', figs);
%
% model_type: 'lumped' | 'discretised' | 'dns' | '2D' | '2DNS'
%   All five share the same 7-station (0..6) TMS field names (P1..P6,
%   P1_0..P6_0, T1..T6, v1..v6, M1..M6, M_dot_2..M_dot_6, m_spill,
%   P_inf, P_inf_tot, T_inf, M_inf, V_inf).
%   'discretised' / 'dns' / '2D' / '2DNS' additionally carry HX-interior
%   arrays (T_air_seg, P_air_seg, P_0_air_seg, v_air_seg, v_channel_seg,
%   m_dot_seg, N_segments), plotted as a smooth curve between stations
%   3 and 4; 'lumped' has no interior resolution and is skipped there.
%
%   For 2D discretised results, most HX-interior arrays are
%   (N_segments+1) x N_seg_cool -- one column per coolant-side segment
%   chain. Only the mean across those columns is plotted, as a single
%   line with a single legend entry (see minmax below for showing the
%   first/last chains too). For N_seg_cool == 1 (old iterative 1D, or 2D
%   trivially run with one column) this is a no-op: the mean of one
%   column is that column. m_dot_seg is the exception -- air mass flow is
%   conserved along the duct, so it's always a plain 1 x (N_segments+1)
%   profile regardless of N_seg_cool, and is plotted as-is (no chains to
%   collapse).
%
% minmax (hard-coded flag, just below): false by default. Set true to
% also plot the first and last coolant-direction segment chains
% (columns 1 and end of each HX-interior array) alongside the mean, as
% thinner unlabeled lines (no extra legend entries) in the same
% interior color/style. Has no visible effect when there's only one
% column (1D, or 2D with N_cool_seg == 1).
%
% Optional name-value pairs:
%   'figs'           - existing figs struct to add to (default: new struct)
%   'color'          - override this call's station marker/line color
%   'interior_color' - override this call's HX-interior line color
%   'marker'         - override this call's station marker
%   'linestyle'      - override this call's station line style
%   'linewidth'      - override this call's LineWidth
%   'markersize'     - override this call's MarkerSize
%   'label'          - override this call's legend label (interior line
%                       is labelled '<label> HX interior')
%
% Sweep pattern, reusing the same figs struct:
%   figs = struct();
%   for i = 1:length(N_list)
%       results_D = run_disc_model_fwdpass(...);
%       figs = plot_stations(results_D, 'discretised', 'figs', figs, 'color', cmap(i,:));
%   end

minmax = false;   % hard-coded; set true to also show first/last coolant-direction chains (unlabeled) for 2D results

p = inputParser;
addParameter(p, 'figs', struct());
addParameter(p, 'color', []);
addParameter(p, 'interior_color', []);
addParameter(p, 'marker', '');
addParameter(p, 'linestyle', '');
addParameter(p, 'linewidth', []);
addParameter(p, 'markersize', []);
addParameter(p, 'label', '');
parse(p, varargin{:});

figs_in   = p.Results.figs;
is_lumped = strcmpi(model_type, 'lumped');

style = get_style(model_type, p.Results.color, p.Results.interior_color, ...
                   p.Results.marker, p.Results.linestyle, p.Results.linewidth, p.Results.markersize);
[lbl, interior_lbl] = get_label(model_type, p.Results.label);

%% ---- Shared freestream quantities ---- %%
P_inf     = results.P_inf;
P_inf_tot = results.P_inf_tot;
T_inf     = results.T_inf;
M_inf     = results.M_inf;
V_inf     = results.V_inf;

%% ---- Station arrays (every model type carries these) ---- %%
P_arr  = [P_inf,     results.P1,   results.P2,   results.P3,   results.P4,   results.P5,   results.P6];
P0_arr = [P_inf_tot, results.P1_0, results.P2_0, results.P3_0, results.P4_0, results.P5_0, results.P6_0];
T_arr  = [T_inf,     results.T1,   results.T2,   results.T3,   results.T4,   results.T5,   results.T6];
v_arr  = [V_inf,     results.v1,   results.v2,   results.v3,   results.v4,   results.v5,   results.v6];
M_arr  = [M_inf,     results.M1,   results.M2,   results.M3,   results.M4,   results.M5,   results.M6];

m_dot_streamtube = results.M_dot_2 + results.m_spill;
M_dot_arr = [m_dot_streamtube, m_dot_streamtube, results.M_dot_2, results.M_dot_3, ...
             results.M_dot_4, results.M_dot_5, results.M_dot_6];

%% ---- HX-interior arrays (discretised-family types only) ---- %%
if ~is_lumped
    x_HX          = linspace(3, 4, results.N_segments+1);
    T_air_seg     = results.T_air_seg;
    P_air_seg     = results.P_air_seg;
    P_0_air_seg   = results.P_0_air_seg;
    v_air_seg     = results.v_air_seg;
    v_channel_seg = results.v_channel_seg;
    m_dot_seg     = results.m_dot_seg;
end

x_stations = 0:6;
labels     = {'0 (\infty)','1','2','3','4','5','6'};

figs_out = figs_in;

    %% ---- Nested: station-level line (marker+line, at the 7 TMS stations) ---- %%
    function add_station_line(ax, y)
        plot(ax, x_stations, y, [style.marker style.linestyle], 'Color', style.color, ...
             'LineWidth', style.lw, 'MarkerFaceColor', style.color, ...
             'MarkerSize', style.ms, 'DisplayName', lbl);
    end

    %% ---- Nested: HX-interior line (smooth curve between stations 3-4) ---- %%
    % Most HX-interior arrays (T_air_seg, P_air_seg, v_air_seg, ...) are
    % (N_segments+1) x N_seg_cool for 2D discretised results -- one column
    % per coolant-direction chain. Others (m_dot_seg: mass flow is
    % conserved along the duct, so it's never split per chain) are always
    % a plain 1 x (N_segments+1) profile, in EITHER discretisation mode.
    % So "multi-chain" is decided by y actually being a matrix, not by
    % which dimension happens to be >1 -- a plain vector is plotted as-is
    % regardless of orientation, and only a true matrix gets the
    % mean/minmax chain-collapsing treatment.
    function add_interior_line(ax, y, ln_style, name_suffix)
        if nargin < 3 || isempty(ln_style), ln_style = '-'; end
        if nargin < 4, name_suffix = ''; end

        if isvector(y)
            y_line = y;
        else
            % Genuine 2D array: orient so length runs along rows (matching
            % x_HX) before collapsing columns, in case it came in transposed.
            if size(y,1) ~= numel(x_HX) && size(y,2) == numel(x_HX)
                y = y.';
            end
            if minmax
                plot(ax, x_HX, y(:,1),   ln_style, 'Color', style.interior_color, ...
                     'LineWidth', max(style.lw*0.5, 0.5), 'HandleVisibility', 'off');
                plot(ax, x_HX, y(:,end), ln_style, 'Color', style.interior_color, ...
                     'LineWidth', max(style.lw*0.5, 0.5), 'HandleVisibility', 'off');
            end
            y_line = mean(y, 2);
        end

        plot(ax, x_HX, y_line, ln_style, 'Color', style.interior_color, 'LineWidth', style.lw, ...
             'DisplayName', [interior_lbl name_suffix]);
    end

%% ---- Figure: Static pressure through TMS ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'static_pressure', 1, 1, 'Static Pressure — TMS Comparison');
if ~isfield(figs_in,'static_pressure')
    set(ax(1), 'XTick', x_stations, 'XTickLabel', labels)
    xlabel(ax(1),'TMS Station'); ylabel(ax(1),'Static pressure [kPa]');
    title(ax(1),'Static pressure through TMS');
end
add_station_line(ax(1), P_arr/1e3);
if ~is_lumped, add_interior_line(ax(1), P_air_seg/1e3); end
legend(ax(1), 'Location','best');
figs_out.static_pressure = fig;

%% ---- Figure: Total pressure through TMS ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'total_pressure', 1, 1, 'Total Pressure — TMS Comparison');
if ~isfield(figs_in,'total_pressure')
    set(ax(1), 'XTick', x_stations, 'XTickLabel', labels)
    xlabel(ax(1),'TMS Station'); ylabel(ax(1),'Total pressure [kPa]');
    title(ax(1),'Total pressure through TMS');
end
add_station_line(ax(1), P0_arr/1e3);
if ~is_lumped, add_interior_line(ax(1), P_0_air_seg/1e3); end
legend(ax(1), 'Location','best');
figs_out.total_pressure = fig;

%% ---- Figure: Static temperature through TMS ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'static_temperature', 1, 1, 'Static Temperature — TMS Comparison');
if ~isfield(figs_in,'static_temperature')
    set(ax(1), 'XTick', x_stations, 'XTickLabel', labels)
    xlabel(ax(1),'TMS Station'); ylabel(ax(1),'Static temperature [K]');
    title(ax(1),'Static temperature through TMS');
end
add_station_line(ax(1), T_arr);
if ~is_lumped, add_interior_line(ax(1), T_air_seg); end
legend(ax(1), 'Location','best');
figs_out.static_temperature = fig;

%% ---- Figure: Velocity through TMS ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'velocity', 1, 1, 'Velocity — TMS Comparison');
if ~isfield(figs_in,'velocity')
    set(ax(1), 'XTick', x_stations, 'XTickLabel', labels)
    xlabel(ax(1),'TMS Station'); ylabel(ax(1),'Velocity [m/s]');
    title(ax(1),'Velocity through TMS');
end
add_station_line(ax(1), v_arr);
if ~is_lumped
    add_interior_line(ax(1), v_air_seg, '-', '');
    add_interior_line(ax(1), v_channel_seg, '--', ' channel');
end
legend(ax(1), 'Location','best');
figs_out.velocity = fig;

%% ---- Figure: Mach number through TMS ---- %%
% No HX-interior resolution tracked for Mach number - station points only.
[fig, ax] = get_or_create_fig(figs_in, 'mach_number', 1, 1, 'Mach Number — TMS Comparison');
if ~isfield(figs_in,'mach_number')
    set(ax(1), 'XTick', x_stations, 'XTickLabel', labels)
    xlabel(ax(1),'TMS Station'); ylabel(ax(1),'Mach number [-]');
    title(ax(1),'Mach number through TMS');
end
add_station_line(ax(1), M_arr);
legend(ax(1), 'Location','best');
figs_out.mach_number = fig;

%% ---- Figure: Combined 2x2 summary panel ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'summary', 2, 2, 'TMS Summary — Model Comparison');
if ~isfield(figs_in,'summary')
    sgtitle(fig, 'TMS Flow Variables — Model Comparison', 'FontSize', 12);
    for i = 1:4
        set(ax(i), 'XTick', x_stations, 'XTickLabel', labels)
    end
    ylabel(ax(1),'Static pressure [kPa]');  title(ax(1),'Static pressure');
    ylabel(ax(2),'Static temperature [K]'); title(ax(2),'Static temperature');
    ylabel(ax(3),'Velocity [m/s]');         title(ax(3),'Velocity');
    ylabel(ax(4),'Mach number [-]');        title(ax(4),'Mach number');
end
add_station_line(ax(1), P_arr/1e3); if ~is_lumped, add_interior_line(ax(1), P_air_seg/1e3); end
add_station_line(ax(2), T_arr);     if ~is_lumped, add_interior_line(ax(2), T_air_seg); end
add_station_line(ax(3), v_arr);     if ~is_lumped, add_interior_line(ax(3), v_air_seg); end
add_station_line(ax(4), M_arr);
for i = 1:4
    legend(ax(i), 'Location','best','FontSize',8);
end
figs_out.summary = fig;

%% ---- Figure: Mass flow through TMS ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'mass_flow_tms', 1, 1, 'Mass Flow — TMS Comparison');
if ~isfield(figs_in,'mass_flow_tms')
    set(ax(1), 'XTick', x_stations, 'XTickLabel', labels)
    xlabel(ax(1),'TMS Station'); ylabel(ax(1),'Mass flow rate [kg/s]');
    title(ax(1),'Mass flow rate through TMS');
end
add_station_line(ax(1), M_dot_arr);
if ~is_lumped, add_interior_line(ax(1), m_dot_seg); end
legend(ax(1), 'Location','best');
figs_out.mass_flow_tms = fig;

end

%% ---- Local helper: default styling per model type ---- %%
function style = get_style(model_type, color_override, interior_color_override, ...
                            marker_override, linestyle_override, lw_override, ms_override)
    switch lower(model_type)
        case 'lumped'
            col = [0.122 0.471 0.706]; icol = [];                      marker = 'o'; ln = '-';  lw = 1.5; ms = 7;
        case 'discretised'
            col = [0.839 0.153 0.157]; icol = [0.173 0.627 0.173];     marker = 's'; ln = '--'; lw = 1.5; ms = 7;
        case 'dns'
            col = [0.580 0.404 0.741]; icol = [1.000 0.498 0.055];     marker = '^'; ln = ':';  lw = 1.5; ms = 7;
        case '2d'
            col = [0.301 0.745 0.933]; icol = [0.494 0.184 0.556];     marker = 'd'; ln = '-.'; lw = 1.5; ms = 7;
        case '2dns'
            col = [0.635 0.078 0.184]; icol = [0.929 0.694 0.125];     marker = 'p'; ln = ':';  lw = 1.5; ms = 7;
        otherwise
            error('plot_stations:UnknownType', ...
                'Unknown model_type "%s". Expected lumped | discretised | dns | 2D | 2DNS.', model_type);
    end

    style.color          = color_override;          if isempty(style.color),          style.color          = col; end
    style.interior_color = interior_color_override;  if isempty(style.interior_color), style.interior_color = icol; end
    style.marker         = marker_override;          if isempty(style.marker),         style.marker         = marker; end
    style.linestyle       = linestyle_override;      if isempty(style.linestyle),      style.linestyle       = ln; end
    style.lw              = lw_override;             if isempty(style.lw),             style.lw              = lw; end
    style.ms              = ms_override;             if isempty(style.ms),             style.ms              = ms; end
end

%% ---- Local helper: default legend labels per model type ---- %%
function [lbl, interior_lbl] = get_label(model_type, label_override)
    if ~isempty(label_override)
        lbl = label_override;
    else
        switch lower(model_type)
            case 'lumped',       lbl = 'Lumped';
            case 'discretised',  lbl = 'Discretised';
            case 'dns',          lbl = 'DNS-based';
            case '2d',           lbl = '2D Discretised';
            case '2dns',         lbl = '2D DNS-based';
        end
    end
    interior_lbl = [lbl ' HX interior'];
end

%% ---- Local helper: get or create a tagged multi-panel figure ---- %%
function [fig, ax] = get_or_create_fig(figs_in, key, nrows, ncols, fig_title)
    if isfield(figs_in, key) && ~isempty(figs_in.(key)) && isvalid(figs_in.(key))
        fig = figs_in.(key);
        figure(fig);
        ax = gobjects(1, nrows*ncols);
        for i = 1:nrows*ncols
            found = findobj(fig, 'Type','axes', 'Tag', sprintf('%s_ax%d', key, i));
            ax(i) = found(1);
        end
    else
        fig = figure('Name', fig_title, 'NumberTitle','off');
        ax = gobjects(1, nrows*ncols);
        for i = 1:nrows*ncols
            ax(i) = subplot(nrows, ncols, i);
            hold(ax(i), 'on');
            grid(ax(i), 'on');
            ax(i).Tag = sprintf('%s_ax%d', key, i);
        end
    end
end