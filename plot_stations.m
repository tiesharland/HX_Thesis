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
%   m_dot_seg, M_air_seg (Mach number), T_0_air_seg, N_segments). The
%   velocity plotted through the HX is v_channel_seg, the actual air
%   velocity in the channels (v_air_seg, the equivalent-duct velocity, is not plotted).
%
%   One line per call: stations 0, 1, 2, then the HX-interior values (which
%   start at station 3 and end at station 4), then stations 5 and 6, all
%   connected. Markers sit on the 7 stations only. 'lumped' has no interior
%   resolution, so it is the plain 7-station line (3 -> 4 straight).
%
%   Default colors match plot_HX (lumped black, discretised red, DNS blue,
%   2D magenta, 2DNS cyan).
%
%   Total temperature at the stations is derived from the stored static
%   temperature and Mach number (T0 = T*(1+(gamma-1)/2*M^2), gamma from
%   gamma_air(cp_air(T))) since the results do not store station totals;
%   the HX interior uses the stored T_0_air_seg.
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
%   'interior_color' - ignored (kept so old calls don't error); each call is
%                       now ONE line in ONE color
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
    v_channel_seg = results.v_channel_seg;
    m_dot_seg     = results.m_dot_seg;
    M_air_seg     = results.M_air_seg;
    T_0_air_seg   = results.T_0_air_seg;
else
    % No interior resolution: empty arrays -> add_station_line draws the
    % plain 7-station line.
    T_air_seg = []; P_air_seg = []; P_0_air_seg = [];
    v_channel_seg = []; m_dot_seg = []; M_air_seg = []; T_0_air_seg = [];
end

% Station total temperature, derived (see header)
gam_arr = arrayfun(@(T) gamma_air(cp_air(T)), T_arr);
T0_arr  = T_arr .* (1 + 0.5*(gam_arr - 1).*M_arr.^2);

x_stations = 0:6;
labels     = {'0 (\infty)','1','2','3','4','5','6'};

figs_out = figs_in;

    %% ---- Nested: one connected line through the TMS ---- %%
    % y_st: the 7 station values (stations 0..6). y_int: HX-interior values
    % (empty for lumped). With an interior, the line runs
    %   station 0, 1, 2, interior (x = 3 ... 4), station 5, 6
    % i.e. the interior replaces the straight 3 -> 4 segment (its first and
    % last points are stations 3 and 4). Markers only on the 7 stations.
    function add_station_line(ax, y_st, y_int)
        if nargin < 3, y_int = []; end
        if isempty(y_int)
            x_line = x_stations;
            y_line = y_st;
            mk_idx = 1:7;
        else
            yi = collapse_chains(y_int);
            n_hx = numel(x_HX);
            x_line = [0 1 2, x_HX(:)', 5 6];
            y_line = [y_st(1:3), yi(:)', y_st(6:7)];
            mk_idx = [1 2 3 4 3+n_hx 4+n_hx 5+n_hx];
            if minmax && ~isvector(y_int)
                add_chain_extremes(ax, y_int, '-');
            end
        end
        plot(ax, x_line, y_line, 'LineStyle', style.linestyle, 'Color', style.color, ...
             'LineWidth', style.lw, 'Marker', style.marker, 'MarkerIndices', mk_idx, ...
             'MarkerFaceColor', style.color, 'MarkerSize', style.ms, 'DisplayName', lbl);
    end

    %% ---- Nested: first/last coolant-direction chain, thin and unlabeled ---- %%
    function add_chain_extremes(ax, y, ln_style)
        if size(y,1) ~= numel(x_HX) && size(y,2) == numel(x_HX), y = y.'; end
        plot(ax, x_HX, y(:,1),   ln_style, 'Color', style.color, ...
             'LineWidth', max(style.lw*0.5, 0.5), 'HandleVisibility', 'off');
        plot(ax, x_HX, y(:,end), ln_style, 'Color', style.color, ...
             'LineWidth', max(style.lw*0.5, 0.5), 'HandleVisibility', 'off');
    end

    %% ---- Nested: collapse a (N+1) x N_seg_cool array to its mean over the coolant direction ---- %%
    % A plain vector (e.g. m_dot_seg, or any 1D-path array) is returned as-is.
    function yc = collapse_chains(y)
        if isvector(y)
            yc = y;
        else
            if size(y,1) ~= numel(x_HX) && size(y,2) == numel(x_HX)
                y = y.';
            end
            yc = mean(y, 2);
        end
    end

%% ---- Figure: Static pressure through TMS ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'static_pressure', 1, 1, 'Static Pressure — TMS Comparison');
if ~isfield(figs_in,'static_pressure')
    set(ax(1), 'XTick', x_stations, 'XTickLabel', labels)
    xlabel(ax(1),'TMS Station'); ylabel(ax(1),'Static pressure [kPa]');
    title(ax(1),'Static pressure through TMS');
end
add_station_line(ax(1), P_arr/1e3, P_air_seg/1e3);
legend(ax(1), 'Location','best');
figs_out.static_pressure = fig;

%% ---- Figure: Total pressure through TMS ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'total_pressure', 1, 1, 'Total Pressure — TMS Comparison');
if ~isfield(figs_in,'total_pressure')
    set(ax(1), 'XTick', x_stations, 'XTickLabel', labels)
    xlabel(ax(1),'TMS Station'); ylabel(ax(1),'Total pressure [kPa]');
    title(ax(1),'Total pressure through TMS');
end
add_station_line(ax(1), P0_arr/1e3, P_0_air_seg/1e3);
legend(ax(1), 'Location','best');
figs_out.total_pressure = fig;

%% ---- Figure: Static temperature through TMS ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'static_temperature', 1, 1, 'Static Temperature — TMS Comparison');
if ~isfield(figs_in,'static_temperature')
    set(ax(1), 'XTick', x_stations, 'XTickLabel', labels)
    xlabel(ax(1),'TMS Station'); ylabel(ax(1),'Static temperature [K]');
    title(ax(1),'Static temperature through TMS');
end
add_station_line(ax(1), T_arr, T_air_seg);
legend(ax(1), 'Location','best');
figs_out.static_temperature = fig;

%% ---- Figure: Total temperature through TMS ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'total_temperature', 1, 1, 'Total Temperature — TMS Comparison');
if ~isfield(figs_in,'total_temperature')
    set(ax(1), 'XTick', x_stations, 'XTickLabel', labels)
    xlabel(ax(1),'TMS Station'); ylabel(ax(1),'Total temperature [K]');
    title(ax(1),'Total temperature through TMS');
end
add_station_line(ax(1), T0_arr, T_0_air_seg);
legend(ax(1), 'Location','best');
figs_out.total_temperature = fig;

%% ---- Figure: Velocity through TMS ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'velocity', 1, 1, 'Velocity — TMS Comparison');
if ~isfield(figs_in,'velocity')
    set(ax(1), 'XTick', x_stations, 'XTickLabel', labels)
    xlabel(ax(1),'TMS Station'); ylabel(ax(1),'Velocity [m/s]');
    title(ax(1),'Velocity through TMS');
end
add_station_line(ax(1), v_arr, v_channel_seg);
legend(ax(1), 'Location','best');
figs_out.velocity = fig;

%% ---- Figure: Mach number through TMS ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'mach_number', 1, 1, 'Mach Number — TMS Comparison');
if ~isfield(figs_in,'mach_number')
    set(ax(1), 'XTick', x_stations, 'XTickLabel', labels)
    xlabel(ax(1),'TMS Station'); ylabel(ax(1),'Mach number [-]');
    title(ax(1),'Mach number through TMS');
end
add_station_line(ax(1), M_arr, M_air_seg);
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
add_station_line(ax(1), P_arr/1e3, P_air_seg/1e3);
add_station_line(ax(2), T_arr,     T_air_seg);
add_station_line(ax(3), v_arr,     v_channel_seg);
add_station_line(ax(4), M_arr,     M_air_seg);
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
add_station_line(ax(1), M_dot_arr, m_dot_seg);
legend(ax(1), 'Location','best');
figs_out.mass_flow_tms = fig;

end

%% ---- Local helper: default styling per model type ---- %%
% Colors are the same as plot_HX's defaults (lumped black, discretised
% autumn(1) = red, DNS winter(1) = blue, 2D spring(1) = magenta, 2DNS
% cool(1) = cyan).
function style = get_style(model_type, color_override, ~, ...
                            marker_override, linestyle_override, lw_override, ms_override)
    switch lower(model_type)
        case 'lumped'
            col = [0 0 0];  marker = 'o'; ln = '-'; lw = 1.5; ms = 7;
        case 'discretised'
            col = autumn(1); col = col(1,:); marker = 's'; ln = '-'; lw = 1.5; ms = 7;
        case 'dns'
            col = winter(1); col = col(1,:); marker = '^'; ln = '-'; lw = 1.5; ms = 7;
        case '2d'
            col = spring(1); col = col(1,:); marker = 'd'; ln = '-'; lw = 1.5; ms = 7;
        case '2dns'
            col = cool(1);   col = col(1,:); marker = 'p'; ln = '-'; lw = 1.5; ms = 7;
        otherwise
            error('plot_stations:UnknownType', ...
                'Unknown model_type "%s". Expected lumped | discretised | dns | 2D | 2DNS.', model_type);
    end

    style.color      = color_override;     if isempty(style.color),      style.color      = col;    end
    style.marker     = marker_override;    if isempty(style.marker),     style.marker     = marker; end
    style.linestyle  = linestyle_override; if isempty(style.linestyle),  style.linestyle  = ln;     end
    style.lw         = lw_override;        if isempty(style.lw),         style.lw         = lw;     end
    style.ms         = ms_override;        if isempty(style.ms),         style.ms         = ms;     end
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