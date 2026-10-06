function figs_out = plot_HX(results, model_type, varargin)
% plot_HX  Add one model's HX results to a shared set of comparison figures.
%
% Usage (one call per dataset - pass the SAME figs struct back in each time
% to keep adding lines to the same figures):
%   figs = plot_HX(results_D,    'discretised');
%   figs = plot_HX(results_L,    'lumped',      'figs', figs);
%   figs = plot_HX(results_DNS,  'dns',         'figs', figs, 'color', col_dns(i,:));
%   figs = plot_HX(results_2D,   '2D',          'figs', figs);
%   figs = plot_HX(results_2DNS, '2DNS',        'figs', figs);
%   figs = plot_HX(results_cfd,  'cfd',         'figs', figs);
%
% model_type: 'lumped' | 'discretised' | 'dns' | '2D' | '2DNS' | 'cfd'
%   - 'lumped' results are scalar (one value per quantity).
%   - 'discretised' / 'dns' / '2D' / '2DNS' all share the same field
%     layout (node arrays sized (N_segments+1) x N_seg_cool, segment
%     arrays sized N_segments x N_seg_cool) and are plotted identically,
%     mean-averaged across the coolant/width direction; they differ only
%     in default styling/label and in the data source that produced them.
%   - 'cfd' is the struct returned by process_wall_CFD_data (fields z, q,
%     Twall, Tbulk, h, Nu, is_entrance, ...). It only has air bulk
%     temperature, h and Nu at discrete z-stations - no pressure, coolant,
%     velocity, Re or Q data - so it only appears on those three subplots
%     (air_state temperature, thermal_perf h, heat_transfer Nu) and is
%     skipped everywhere else. Position is plotted as -z (CFD flow runs in
%     -z; the model's x runs in the flow direction from 0), in mm, and the
%     z=0 leading-edge station is excluded (h formally diverges there for
%     a step-change wall condition, so it's a mesh-dependent artifact, not
%     a converged value - see process_wall_CFD_data's own exclusion of it
%     from summary averages).
%
% Optional name-value pairs:
%   'figs'       - existing figs struct to add to (default: new struct)
%   'color'      - override this call's line/marker color
%   'marker'     - override this call's marker (e.g. 'o', '^', 's')
%   'linestyle'  - override this call's line style (e.g. '--', ':', '-.')
%   'linewidth'  - override this call's LineWidth
%   'markersize' - override this call's MarkerSize
%   'label'      - override this call's legend label
%
% Sweep pattern (e.g. over N_segments), reusing the same figs struct:
%   figs = struct();
%   cmap = autumn(length(N_list));
%   for i = 1:length(N_list)
%       results_D = run_disc_model_fwdpass(...);
%       figs = plot_HX(results_D, 'discretised', 'figs', figs, 'color', cmap(i,:));
%   end

p = inputParser;
addParameter(p, 'figs', struct());
addParameter(p, 'color', []);
addParameter(p, 'marker', '');
addParameter(p, 'linestyle', '');
addParameter(p, 'linewidth', []);
addParameter(p, 'markersize', []);
addParameter(p, 'label', '');
parse(p, varargin{:});

figs_in = p.Results.figs;

is_lumped = strcmpi(model_type, 'lumped');
is_cfd    = strcmpi(model_type, 'cfd');
is_family = ~is_lumped && ~is_cfd;   % discretised / dns / 2D / 2DNS

style = get_style(model_type, p.Results.color, p.Results.marker, p.Results.linewidth, p.Results.markersize, p.Results.linestyle);
lbl   = get_label(results, model_type, p.Results.label);

% Shared x-axis positions [mm]
if is_lumped
    xb = [0, results.L_solution] * 1000;   % two-point flat line, no midpoints
    xm = [];
elseif is_cfd
    cfd_valid = ~results.is_entrance;      % exclude the z=0 entrance singularity
    x_cfd = -results.z(cfd_valid) * 1000;  % flip sign (CFD flow runs in -z), m -> mm
    xb = []; xm = [];                      % unused for cfd
else
    xb = linspace(0, results.L_solution, results.N_segments+1) * 1000;
    xm = 0.5*(xb(1:end-1) + xb(2:end));
end

% Not every model_type's results struct carries T_h_i/T_c_i (the CFD
% struct doesn't) - guard the reference-line setup below accordingly.
if isfield(results, 'T_h_i'), T_h_i = results.T_h_i; else, T_h_i = []; end
if isfield(results, 'T_c_i'), T_c_i = results.T_c_i; else, T_c_i = []; end

figs_out = figs_in;

    %% ---- Nested: add one styled line to an axes, using this call's style/label ---- %%
    function ln = addline(ax, x, y, varargin_extra)
        if nargin < 4, varargin_extra = {}; end
        ln = plot(ax, x, y, [style.line style.marker], 'Color', style.color, ...
                  'LineWidth', style.lw, 'MarkerFaceColor', style.color, ...
                  'MarkerSize', style.ms, 'DisplayName', lbl, varargin_extra{:});
    end

%% ---- Figure: Air State Variables ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'air_state', 2, 1, 'HX Air State Variables');

if ~isfield(figs_in,'air_state')
    if ~isempty(T_h_i), yline(ax(1), T_h_i, '--r', 'DisplayName', 'T_{h,i}', 'LabelHorizontalAlignment','left'); end
    if ~isempty(T_c_i), yline(ax(1), T_c_i, '--b', 'DisplayName', 'T_{c,i}', 'LabelHorizontalAlignment','left'); end
    xlabel(ax(1),'Position along HX [mm]'); ylabel(ax(1),'Air temperature [K]');
    title(ax(1),'Air temperature distribution through HX');
end
if is_lumped
    addline(ax(1), xb, [results.T3, results.T4]);
    addline(ax(1), xb, [results.T_mean_c, results.T_mean_c], {'LineStyle','--', 'DisplayName','Lumped T_{mean,c}'});
elseif is_cfd
    addline(ax(1), x_cfd, results.Tbulk(cfd_valid));   % CFD air bulk temperature
elseif is_family
    addline(ax(1), xb, mean(results.T_air_seg,2));
end
legend(ax(1), 'Location','southeast');

if ~isfield(figs_in,'air_state')
    xlabel(ax(2),'Position along HX [mm]'); ylabel(ax(2),'Air static pressure [kPa]');
    title(ax(2),'Air pressure distribution through HX');
end
if is_lumped
    addline(ax(2), xb, [results.P3, results.P4]/1000);
elseif is_family
    addline(ax(2), xb, mean(results.P_air_seg,2)/1000);
end
legend(ax(2), 'Location','northeast');
figs_out.air_state = fig;

%% ---- Figure: Coolant Variables ---- %%
% No coolant data in the CFD struct - skipped entirely for 'cfd'.
[fig, ax] = get_or_create_fig(figs_in, 'coolant', 3, 1, 'HX Coolant Variables');

if ~isfield(figs_in,'coolant')
    if ~isempty(T_h_i), yline(ax(1), T_h_i, '--r', 'DisplayName', 'T_{h,i}', 'LabelHorizontalAlignment','left'); end
    if ~isempty(T_c_i), yline(ax(1), T_c_i, '--b', 'DisplayName', 'T_{c,i}', 'LabelHorizontalAlignment','left'); end
    xlabel(ax(1),'Position along HX [mm]'); ylabel(ax(1),'Coolant temperature [K]');
    title(ax(1),'Coolant outlet temperature distribution through HX');
end
if is_lumped
    addline(ax(1), xb, [results.T_cool_out, results.T_cool_out]);
elseif is_family
    addline(ax(1), xm, results.T_cool_seg(:,end));
end
legend(ax(1), 'Location','southeast');

if ~isfield(figs_in,'coolant')
    xlabel(ax(2),'Position along HX [mm]'); ylabel(ax(2),'Coolant pressure drop [bar]');
    title(ax(2),'Coolant pressure drop distribution through HX');
end
if is_lumped
    addline(ax(2), xb, [results.dp_coolant, results.dp_coolant]/1e5);
elseif is_family
    addline(ax(2), xm, mean(results.dp_cool_seg,2)/1e5);
end
legend(ax(2), 'Location','southeast');

if ~isfield(figs_in,'coolant')
    xlabel(ax(3),'Position along HX [mm]'); ylabel(ax(3),'Temperature difference [K]');
    title(ax(3),'(Mean) Temperature difference distribution through HX');
end
if is_lumped
    dT = results.T_mean_h - results.T_mean_c;
    addline(ax(3), xb, [dT, dT]);
elseif is_family
    dT_seg = mean(movmean(results.T_cool_seg,2,2,'Endpoints','discard') - ...
                  movmean(results.T_air_seg,2,'Endpoints','discard'), 2);
    addline(ax(3), xm, dT_seg);
end
legend(ax(3), 'Location','southeast');
figs_out.coolant = fig;

%% ---- Figure: Coolant temperature across the width (first vs last air-side segment) ---- %%
% T_cool_seg is N_segments x (N_seg_cool+1): column 1 is the coolant inlet
% temperature, column l+1 the coolant temperature after width-wise cell l.
% Top: first air-side segment; bottom: last.
%   2D path (N_cool_seg >= 1): just the coolant temperature across the
%       width. (A cell's mean temperature T_mean_h_arr(k,l) is simply the
%       average of its inlet/outlet nodes, so it is not drawn separately.)
%   1D path (N_cool_seg == 0): the coolant has no width resolution, so the
%       curve is just the two nodes (inlet, outlet) of that segment.
if is_family
    [fig, ax] = get_or_create_fig(figs_in, 'coolant_width', 2, 1, 'Coolant Temperature Across Width');

    n_air_seg = size(results.T_cool_seg, 1);
    x_w       = linspace(0, 1, size(results.T_cool_seg, 2));
    seg_idx   = [1, n_air_seg];
    seg_names = {'First air-side segment', 'Last air-side segment'};

    for j = 1:2
        if ~isfield(figs_in,'coolant_width')
            if ~isempty(T_h_i), yline(ax(j), T_h_i, '--r', 'DisplayName', 'T_{h,i}', 'LabelHorizontalAlignment','left'); end
            xlabel(ax(j),'Position across coolant width [-]  (0 = coolant inlet)'); ylabel(ax(j),'Coolant temperature [K]');
            title(ax(j), seg_names{j});
        end
        addline(ax(j), x_w, results.T_cool_seg(seg_idx(j), :));
        legend(ax(j), 'Location','best');
    end
    figs_out.coolant_width = fig;
end

%% ---- Figure: Coolant-side thermal performance (h, q) ---- %%
% Same calculation as the air-side 'HX Thermal Performance' figure, vs
% position along the HX (width-averaged per air-side segment):
%   h : coolant-side heat transfer coefficient (results.h_cool_seg)
%   q : h_cool*(T_cool - T_air), cell-averaged, same as the air-side q
% Only drawn if the results contain h_cool_seg (results saved before it was
% stored don't have it).
if is_family && isfield(results, 'h_cool_seg')
    [fig, ax] = get_or_create_fig(figs_in, 'coolant_thermal', 2, 1, 'HX Coolant Thermal Performance');
    if ~isfield(figs_in,'coolant_thermal')
        xlabel(ax(1),'Position along HX [mm]'); ylabel(ax(1),'h [W/m^2K]');
        title(ax(1),'Coolant-side heat transfer coefficient through HX');
        xlabel(ax(2),'Position along HX [mm]'); ylabel(ax(2),'q [kW/m^2]');
        title(ax(2),'Coolant-side heat flux through HX');
    end
    addline(ax(1), xm, mean(results.h_cool_seg, 2));
    q_cool = mean(results.h_cool_seg.*(movmean(results.T_cool_seg,2,2,'Endpoints','discard') - ...
                  movmean(results.T_air_seg,2,'Endpoints','discard')), 2)/1000;
    addline(ax(2), xm, q_cool);
    for j = 1:2, legend(ax(j), 'Location','best'); end
    figs_out.coolant_thermal = fig;
end

%% ---- Figure: Cumulative heat rejection across the coolant width (first vs last air-side segment) ---- %%
% Q_seg_arr(k,l) is the heat transferred in width-wise cell l of air-side
% segment k. Plotted cumulatively across the width, like the air-side
% cumulative heat rejection along the length: 0 at the coolant inlet side,
% the segment's total heat transfer at the far side.
%   2D (N_cool_seg >= 1): one point per cell edge (cumulative sum of cells).
%   1D (N_cool_seg == 0): a single cell, so a straight line from 0 to the
%       segment's total Q.
% Top: first air-side segment; bottom: last.
if is_family
    [fig, ax] = get_or_create_fig(figs_in, 'q_width', 2, 1, 'Cumulative Heat Rejection Across Width');

    n_air_seg = size(results.Q_seg_arr, 1);
    n_cells   = size(results.Q_seg_arr, 2);
    seg_idx   = [1, n_air_seg];
    seg_names = {'First air-side segment', 'Last air-side segment'};
    edges     = linspace(0, 1, n_cells+1);

    for j = 1:2
        if ~isfield(figs_in,'q_width')
            xlabel(ax(j),'Position across coolant width [-]  (0 = coolant inlet)');
            ylabel(ax(j),'\Delta Q [kW]');
            title(ax(j), seg_names{j});
        end
        addline(ax(j), edges, [0, cumsum(results.Q_seg_arr(seg_idx(j), :))]/1000);
        legend(ax(j), 'Location','best');
    end
    figs_out.q_width = fig;
end

%% ---- Figure: Air Flow Variables (v, Re, Pr) ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'flow', 3, 1, 'HX Air Flow Variables');

if ~isfield(figs_in,'flow')
    xlabel(ax(1),'Position along HX [mm]'); ylabel(ax(1),'Air channel velocity [m/s]');
    title(ax(1),'Channel air velocity distribution through HX');
end
if is_lumped
    addline(ax(1), xb, [results.v_channel_air, results.v_channel_air]);
elseif is_family
    addline(ax(1), xb, mean(results.v_channel_seg,2));
end
legend(ax(1), 'Location','northeast');

if ~isfield(figs_in,'flow')
    yline(ax(2), 2300, '--k', 'DisplayName', 'Re = 2300', 'LabelHorizontalAlignment','left');
    yline(ax(2), 3000, '--r', 'DisplayName', 'Re = 3000', 'LabelHorizontalAlignment','left');
    xlabel(ax(2),'Position along HX [mm]'); ylabel(ax(2),'Reynolds number [-]');
    title(ax(2),'Air Reynolds number through HX');
end
if is_lumped
    addline(ax(2), xb, [results.Re_air, results.Re_air]);
elseif is_family
    addline(ax(2), xb, mean(results.Re_air_seg,2));
end
legend(ax(2), 'Location','northeast');

if ~isfield(figs_in,'flow')
    xlabel(ax(3),'Position along HX [mm]'); ylabel(ax(3),'Prandtl number [-]');
    title(ax(3),'Air Prandtl number through HX');
end
if is_lumped
    addline(ax(3), xb, [results.Pr_air, results.Pr_air]);
elseif is_family
    addline(ax(3), xb, mean(results.Pr_air_seg,2));  % NOTE: fixed missing dim-2 arg present in the old DNS branch
end
legend(ax(3), 'Location','northeast');
figs_out.flow = fig;

%% ---- Figure: Heat Transfer Variables (f_air, f_hx, Nu) ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'heat_transfer', 3, 1, 'HX Heat Transfer Variables');

if ~isfield(figs_in,'heat_transfer')
    xlabel(ax(1),'Position along HX [mm]'); ylabel(ax(1),'Friction factor f [-]');
    title(ax(1),'Air-side Konakov friction factor through HX');
end
if is_lumped
    addline(ax(1), xb, [results.f_air, results.f_air]);
elseif is_family
    addline(ax(1), xb, mean(results.f_air_seg,2));
end
legend(ax(1), 'Location','southeast');

if ~isfield(figs_in,'heat_transfer')
    xlabel(ax(2),'Position along HX [mm]'); ylabel(ax(2),'Friction factor f [-]');
    title(ax(2),'Air-side Zanke friction factor through HX');
end
if is_lumped
    addline(ax(2), xb, [results.f_hx, results.f_hx]);
elseif is_family
    addline(ax(2), xb, mean(results.f_hx_seg,2));  % NOTE: previously not plotted for DNS - now included for every family type
end
legend(ax(2), 'Location','southeast');

if ~isfield(figs_in,'heat_transfer')
    xlabel(ax(3),'Position along HX [mm]'); ylabel(ax(3),'Nusselt number [-]');
    title(ax(3),'Nusselt number through HX');
end
if is_lumped
    addline(ax(3), xb, [results.Nu_air, results.Nu_air]);
elseif is_cfd
    addline(ax(3), x_cfd, results.Nu(cfd_valid));   % CFD Nusselt number
elseif is_family
    addline(ax(3), xm, mean(results.Nu_air_seg,2));
end
legend(ax(3), 'Location','northeast');
figs_out.heat_transfer = fig;

%% ---- Figure: Thermal Performance (h, q, cumulative Q) ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'thermal_perf', 3, 1, 'HX Thermal Performance');

if ~isfield(figs_in,'thermal_perf')
    xlabel(ax(1),'Position along HX [mm]'); ylabel(ax(1),'h [W/m^2K]');
    title(ax(1),'Air-side heat transfer coefficient through HX');
end
if is_lumped
    addline(ax(1), xb, [results.h_air, results.h_air]);
elseif is_cfd
    addline(ax(1), x_cfd, results.h(cfd_valid));   % CFD heat transfer coefficient
elseif is_family
    addline(ax(1), xm, mean(results.h_air_seg,2));
end
legend(ax(1), 'Location','northeast');

if ~isfield(figs_in,'thermal_perf')
    xlabel(ax(2),'Position along HX [mm]'); ylabel(ax(2),'q [kW/m^2]');
    title(ax(2),'Air-side heat flux through HX');
end
if is_lumped
    q = results.h_air*(mean([T_h_i, results.T_cool_out]) - results.T_mean_c)/1000;
    addline(ax(2), xb, [q, q]);
elseif is_family
    q_seg = mean(results.h_air_seg.*(movmean(results.T_cool_seg,2,2,'Endpoints','discard') - ...
                 movmean(results.T_air_seg,2,'Endpoints','discard')), 2)/1000;
    addline(ax(2), xm, q_seg);
end
legend(ax(2), 'Location','northeast');

if ~isfield(figs_in,'thermal_perf')
    xlabel(ax(3),'Position along HX [mm]'); ylabel(ax(3),'\Delta Q [kW]');
    title(ax(3),'Cumulative heat rejection through HX');
end
if is_lumped
    addline(ax(3), xb, [0, results.Q_pred_solution]/1000);
elseif is_family
    addline(ax(3), xm, cumsum(sum(results.Q_seg_arr,2))/1000);
end
legend(ax(3), 'Location','northeast');
figs_out.thermal_perf = fig;

%% ---- Figure: Pressure Drop ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'pressure_drop', 2, 1, 'HX Pressure Drop');

if ~isfield(figs_in,'pressure_drop')
    xlabel(ax(1),'Segment midpoint position [mm]'); ylabel(ax(1),'\Delta p per segment [Pa]');
    title(ax(1),'Pressure drop per segment');
end
if is_family
    % No per-segment breakdown exists for lumped or cfd; skip ax(1) for those.
    addline(ax(1), xm, mean(results.dp_seg_arr,2));
    legend(ax(1), 'Location','northeast');
end

if ~isfield(figs_in,'pressure_drop')
    xlabel(ax(2),'Position along HX length [mm]'); ylabel(ax(2),'Cumulative \Delta p [Pa]');
    title(ax(2),'Cumulative pressure drop through HX');
end
if is_lumped
    addline(ax(2), xb, [0, results.dp_hx]);
elseif is_family
    addline(ax(2), xb, [0; cumsum(mean(results.dp_seg_arr,2))]);
end
legend(ax(2), 'Location','northwest');
figs_out.pressure_drop = fig;

%% ---- Figure: Mass Flow Conservation ---- %%
% No lumped or cfd equivalent - mass flow is conserved by construction in
% the lumped model, and the CFD struct doesn't carry a mass-flow field.
if is_family
    [fig, ax] = get_or_create_fig(figs_in, 'mass_flow', 1, 1, 'Mass Flow Conservation');

    if ~isfield(figs_in,'mass_flow')
        xlabel(ax(1),'Position along HX [mm]'); ylabel(ax(1),'Air mass flow rate [kg/s]');
        title(ax(1),'Air mass flow rate conservation through HX');
    end
    addline(ax(1), xb, results.m_dot_seg);
    legend(ax(1), 'Location','best');
    figs_out.mass_flow = fig;

    fprintf("\n=== MASS FLOW CONSERVATION (%s, N=%d) ===\n", model_type, results.N_segments);
    fprintf("M_dot_3             = %.6f kg/s\n", results.M_dot_3);
    fprintf("M_dot_4             = %.6f kg/s\n", results.M_dot_4);
    fprintf("Max deviation       = %.2e kg/s\n", max(abs(results.m_dot_seg - results.M_dot_3)));
end

%% ---- Figure: Boundary layer evolution ---- %%
% No lumped or cfd equivalent - boundary layer growth is only tracked
% per segment in the family models.
if is_family
    [fig, ax] = get_or_create_fig(figs_in, 'boundary_layer', 2, 1, 'Boundary Layer Analysis');

    if ~isfield(figs_in,'boundary_layer')
        xlabel(ax(1),'Segment midpoint position [mm]'); ylabel(ax(1),'Boundary layer height [mm]');
        title(ax(1),'Boundary layer height through HX');
    end
    addline(ax(1), xm, mean(results.delta_BL_seg,2)*1000);
    legend(ax(1), 'Location','northeast');

    if ~isfield(figs_in,'boundary_layer')
        xlabel(ax(2),'Segment midpoint position [mm]'); ylabel(ax(2),'Hydraulic diameter [mm]');
        title(ax(2),'Hydraulic diameter through HX');
    end
    addline(ax(2), xm, mean(results.d_h_bulk_seg,2)*1000);
    legend(ax(2), 'Location','northwest');
    figs_out.boundary_layer = fig;   % FIX: previously assigned to figs_out.pressure_drop by mistake
end

end

%% ---- Local helper: default styling per model type ---- %%
function style = get_style(model_type, color_override, marker_override, lw_override, ms_override, ls_override)
    switch lower(model_type)
        case 'lumped'
            cmap = [0 0 0]; marker = 's'; lw = 1.5; ms = 6; ln = '-';
        case 'discretised'
            cmap = autumn(1); cmap = cmap(1,:); marker = 'o'; lw = 0.9; ms = 3; ln = '-';
        case 'dns'
            cmap = winter(1); cmap = cmap(1,:); marker = '^'; lw = 0.9; ms = 3; ln = '-';
        case '2d'
            cmap = spring(1); cmap = cmap(1,:); marker = 'd'; lw = 0.9; ms = 3; ln = '-';
        case '2dns'
            cmap = cool(1); cmap = cmap(1,:); marker = 'p'; lw = 0.9; ms = 3; ln = '-';
        case 'cfd'
            cmap = [0.3 0.3 0.3]; marker = 'x'; lw = 1.0; ms = 5; ln = '-';
        otherwise
            error('plot_HX:UnknownType', ...
                'Unknown model_type "%s". Expected lumped | discretised | dns | 2D | 2DNS | cfd.', model_type);
    end

    style.color  = color_override;  if isempty(style.color),  style.color  = cmap;   end
    style.marker = marker_override; if isempty(style.marker), style.marker = marker; end
    style.lw     = lw_override;     if isempty(style.lw),     style.lw     = lw;     end
    style.ms     = ms_override;     if isempty(style.ms),     style.ms     = ms;     end
    style.line   = ln;  if ~isempty(ls_override), style.line = ls_override; end
end

%% ---- Local helper: default legend label per model type ---- %%
function lbl = get_label(results, model_type, label_override)
    if ~isempty(label_override)
        lbl = label_override;
        return
    end
    switch lower(model_type)
        case 'lumped'
            lbl = sprintf('Lumped (L=%.3fm)', results.L_solution);
        case 'discretised'
            lbl = sprintf('Discretised (N=%d, L=%.3fm)', results.N_segments, results.L_solution);
        case 'dns'
            lbl = sprintf('DNS-based (N=%d, L=%.3fm)', results.N_segments, results.L_solution);
        case '2d'
            if isfield(results, 'N_seg_cool')
                lbl = sprintf('2D Discretised (N=%dx%d, L=%.3fm)', results.N_segments, results.N_seg_cool, results.L_solution);
            else
                lbl = sprintf('2D Discretised (N=%d, L=%.3fm)', results.N_segments, results.L_solution);
            end
        case '2dns'
            if isfield(results, 'N_seg_cool')
                lbl = sprintf('2D DNS-based (N=%dx%d, L=%.3fm)', results.N_segments, results.N_seg_cool, results.L_solution);
            else
                lbl = sprintf('2D DNS-based (N=%d, L=%.3fm)', results.N_segments, results.L_solution);
            end
        case 'cfd'
            lbl = 'CFD (unit cell)';
    end
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