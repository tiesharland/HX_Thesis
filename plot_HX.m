function figs_out = plot_HX(results_D, varargin)
% Usage:
%   figs = plot_HX(results_D);
%   figs = plot_HX(results_D2, 'figs', figs, 'color', c2);
%   figs = plot_HX(results_D, 'lumped', results_L, 'figs', figs, 'plot_lumped', true);
%
% Suggested N_segments sweep call pattern (discretised):
%   colors = summer(length(N_list));
%   for i = 1:length(N_list)
%       figs = plot_HX(results_D_arr(i), 'lumped', results_L, 'figs', figs, ...
%                       'color', colors(i,:), 'plot_lumped', (i==1));
%   end
%
% Suggested N_segments sweep call pattern (DNS, if explored later):
%   colors_dns = cool(length(N_list_DNS));
%   for i = 1:length(N_list_DNS)
%       figs = plot_HX(results_D, 'lumped', results_L, 'dns', results_DNS_arr(i), ...
%                       'figs', figs, 'color_dns', colors_dns(i,:), ...
%                       'plot_lumped', (i==1), 'plot_dns', true);
%   end

p = inputParser;
addParameter(p, 'lumped', []);
addParameter(p, 'dns', []);
addParameter(p, 'figs', struct());
addParameter(p, 'color', []);      % discretised color (per-call, e.g. summer colormap)
addParameter(p, 'color_dns', []);  % dns color (per-call, e.g. cool colormap)
addParameter(p, 'color_lump', [])
addParameter(p, 'label', '');
addParameter(p, 'label_dns', '');
addParameter(p, 'plot_lumped', true);
addParameter(p, 'plot_dns', true);
parse(p, varargin{:});

has_lumped = ~isempty(p.Results.lumped) && p.Results.plot_lumped;
has_dns    = ~isempty(p.Results.dns) && p.Results.plot_dns;
figs_in    = p.Results.figs;

if isempty(p.Results.color)
    sm = autumn(1);
    col_disc = sm(1,:);
else
    col_disc = p.Results.color;
end
lw_disc = 0.9;
ms_disc = 3;

% Lumped: fixed styling/color, only ever drawn once
lw = 1.5; ms = 6;
% col_lump = [0.122 0.471 0.706];
if isempty(p.Results.color_lump)
    % col_lump = [0.122 0.471 0.706];
    col_lump = [0 0 0];
else
    col_lump = p.Results.color_lump;
end 

% DNS: colormap-driven, thin/small, sweepable (mirrors discretised treatment)
if isempty(p.Results.color_dns)
    pm = winter(1);
    col_dns = pm(1,:);
else
    col_dns = p.Results.color_dns;
end
lw_dns = 0.9;
ms_dns = 3;

if has_lumped, results_L = p.Results.lumped; end
if has_dns,    results_DNS = p.Results.dns;  end

%% ---- Only keep locals that are reused across multiple plots/labels ---- %%
N_segments   = results_D.N_segments;
L_solution_D = results_D.L_solution;
T_h_i        = results_D.T_h_i;
T_c_i        = results_D.T_c_i;
M_dot_3      = results_D.M_dot_3;
M_dot_4      = results_D.M_dot_4;

if isempty(p.Results.label)
    lbl = sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D);
else
    lbl = p.Results.label;
end
if has_lumped
    lbl_L = sprintf('Lumped (L=%.3fm)', results_L.L_solution);
end
if has_dns
    if isempty(p.Results.label_dns)
        lbl_DNS = sprintf('DNS-based (N=%d, L=%.3fm)', results_DNS.N_segments, results_DNS.L_solution);
    else
        lbl_DNS = p.Results.label_dns;
    end
end

%% ---- Axis positions ---- %%
x_bounds_D = linspace(0, L_solution_D, N_segments+1);
x_mid_D    = 0.5*(x_bounds_D(1:end-1) + x_bounds_D(2:end));
if has_lumped
    x_bounds_L = [0, results_L.L_solution];
end
if has_dns
    x_bounds_DNS = linspace(0, results_DNS.L_solution, results_DNS.N_segments+1);
    x_mid_DNS    = 0.5*(x_bounds_DNS(1:end-1) + x_bounds_DNS(2:end));
end

figs_out = figs_in;

%% ---- Figure: Air State Variables ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'air_state', 2, 1, 'HX Air State Variables');

if ~isfield(figs_in,'air_state')
    yline(ax(1), T_h_i, '--r', 'DisplayName', 'T_{h,i}', 'LabelHorizontalAlignment','left');
    yline(ax(1), T_c_i, '--b', 'DisplayName', 'T_{c,i}', 'LabelHorizontalAlignment','left');
    xlabel(ax(1),'Position along HX [mm]'); ylabel(ax(1),'Air temperature [K]');
    title(ax(1),'Air temperature distribution through HX');
end
if has_lumped
    plot(ax(1), x_bounds_L*1000, [results_L.T3, results_L.T4], 's-', 'Color', col_lump, 'LineWidth', lw, ...
         'MarkerFaceColor', col_lump, 'MarkerSize', ms, 'DisplayName', lbl_L);
    plot(ax(1), x_bounds_L*1000, [results_L.T_mean_c, results_L.T_mean_c], 's--', 'Color', col_lump, ...
         'LineWidth', lw, 'DisplayName', 'Lumped T_{mean,c}');
end
if has_dns
    plot(ax(1), x_bounds_DNS*1000, results_DNS.T_air_seg, '^-', 'Color', col_dns, 'LineWidth', lw_dns, ...
         'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, 'DisplayName', lbl_DNS);
end
plot(ax(1), x_bounds_D*1000, results_D.T_air_seg, 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(1), 'Location','southeast');

if ~isfield(figs_in,'air_state')
    xlabel(ax(2),'Position along HX [mm]'); ylabel(ax(2),'Air static pressure [kPa]');
    title(ax(2),'Air pressure distribution through HX');
end
if has_lumped
    plot(ax(2), x_bounds_L*1000, [results_L.P3, results_L.P4]/1000, 's-', 'Color', col_lump, 'LineWidth', lw, ...
         'MarkerFaceColor', col_lump, 'MarkerSize', ms, 'DisplayName', lbl_L);
end
if has_dns
    plot(ax(2), x_bounds_DNS*1000, results_DNS.P_air_seg/1000, '^-', 'Color', col_dns, 'LineWidth', lw_dns, ...
         'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, 'DisplayName', lbl_DNS);
end
plot(ax(2), x_bounds_D*1000, results_D.P_air_seg/1000, 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(2), 'Location','northeast');
figs_out.air_state = fig;

%% ---- Figure: Coolant Variables ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'coolant', 2, 1, 'HX Coolant Variables');

if ~isfield(figs_in,'coolant')
    yline(ax(1), T_h_i, '--r', 'DisplayName', 'T_{h,i}', 'LabelHorizontalAlignment','left');
    yline(ax(1), T_c_i, '--b', 'DisplayName', 'T_{c,i}', 'LabelHorizontalAlignment','left');
    xlabel(ax(1),'Position along HX [mm]'); ylabel(ax(1),'Coolant temperature [K]');
    title(ax(1),'Coolant temperature distribution through HX');
end
if has_lumped
    plot(ax(1), x_bounds_L*1000, [results_L.T_cool_out, results_L.T_cool_out], 's-', 'Color', col_lump, ...
         'LineWidth', lw, 'MarkerFaceColor', col_lump, 'MarkerSize', ms, 'DisplayName', lbl_L);
end
if has_dns
    plot(ax(1), x_mid_DNS*1000, results_DNS.T_cool_out_arr, '^-', 'Color', col_dns, 'LineWidth', lw_dns, ...
         'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, 'DisplayName', lbl_DNS);
end
plot(ax(1), x_mid_D*1000, results_D.T_cool_out_arr, 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(1), 'Location','southeast');

if ~isfield(figs_in,'coolant')
    xlabel(ax(2),'Position along HX [mm]'); ylabel(ax(2),'Coolant pressure drop [bar]');
    title(ax(2),'Coolant pressure drop distribution through HX');
end
if has_lumped
    plot(ax(2), x_bounds_L*1000, [results_L.dp_coolant, results_L.dp_coolant]/1e5, 's-', 'Color', col_lump, ...
         'LineWidth', lw, 'MarkerFaceColor', col_lump, 'MarkerSize', ms, 'DisplayName', lbl_L);
end
if has_dns
    plot(ax(2), x_mid_DNS*1000, results_DNS.dp_cool_seg/1e5, '^-', 'Color', col_dns, 'LineWidth', lw_dns, ...
         'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, 'DisplayName', lbl_DNS);
end
plot(ax(2), x_mid_D*1000, results_D.dp_cool_seg/1e5, 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(2), 'Location','southeast');
figs_out.coolant = fig;

%% ---- Figure: Air Flow Variables (v, Re, Pr) ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'flow', 3, 1, 'HX Air Flow Variables');

if ~isfield(figs_in,'flow')
    xlabel(ax(1),'Position along HX [mm]'); ylabel(ax(1),'Air channel velocity [m/s]');
    title(ax(1),'Channel air velocity distribution through HX');
end
if has_lumped
    plot(ax(1), x_bounds_L*1000, [results_L.v_channel_air, results_L.v_channel_air], 's-', 'Color', col_lump, ...
         'LineWidth', lw, 'MarkerFaceColor', col_lump, 'MarkerSize', ms, 'DisplayName', lbl_L);
end
if has_dns
    plot(ax(1), x_bounds_DNS*1000, results_DNS.v_channel_seg, '^-', 'Color', col_dns, 'LineWidth', lw_dns, ...
         'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, 'DisplayName', lbl_DNS);
end
plot(ax(1), x_bounds_D*1000, results_D.v_channel_seg, 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(1), 'Location','northeast');

if ~isfield(figs_in,'flow')
    yline(ax(2), 2300, '--k', 'DisplayName', 'Re = 2300', 'LabelHorizontalAlignment','left');
    yline(ax(2), 3000, '--r', 'DisplayName', 'Re = 3000', 'LabelHorizontalAlignment','left');
    xlabel(ax(2),'Position along HX [mm]'); ylabel(ax(2),'Reynolds number [-]');
    title(ax(2),'Air Reynolds number through HX');
end
if has_lumped
    plot(ax(2), x_bounds_L*1000, [results_L.Re_air, results_L.Re_air], 's-', 'Color', col_lump, 'LineWidth', lw, ...
         'MarkerFaceColor', col_lump, 'MarkerSize', ms, 'DisplayName', lbl_L);
end
if has_dns
    plot(ax(2), x_bounds_DNS*1000, results_DNS.Re_air_seg, '^-', 'Color', col_dns, 'LineWidth', lw_dns, ...
         'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, 'DisplayName', lbl_DNS);
end
plot(ax(2), x_bounds_D*1000, results_D.Re_air_seg, 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(2), 'Location','northeast');

if ~isfield(figs_in,'flow')
    xlabel(ax(3),'Position along HX [mm]'); ylabel(ax(3),'Prandtl number [-]');
    title(ax(3),'Air Prandtl number through HX');
end
if has_lumped
    plot(ax(3), x_bounds_L*1000, [results_L.Pr_air, results_L.Pr_air], 's-', 'Color', col_lump, 'LineWidth', lw, ...
         'MarkerFaceColor', col_lump, 'MarkerSize', ms, 'DisplayName', lbl_L);
end
if has_dns
    plot(ax(3), x_bounds_DNS*1000, results_DNS.Pr_air_seg, '^-', 'Color', col_dns, 'LineWidth', lw_dns, ...
         'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, 'DisplayName', lbl_DNS);
end
plot(ax(3), x_bounds_D*1000, results_D.Pr_air_seg, 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(3), 'Location','northeast');
figs_out.flow = fig;

%% ---- Figure: Heat Transfer Variables (f_air, f_hx, Nu) ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'heat_transfer', 3, 1, 'HX Heat Transfer Variables');

if ~isfield(figs_in,'heat_transfer')
    xlabel(ax(1),'Position along HX [mm]'); ylabel(ax(1),'Friction factor f [-]');
    title(ax(1),'Air-side Konakov friction factor through HX');
end
if has_lumped
    plot(ax(1), x_bounds_L*1000, [results_L.f_air, results_L.f_air], 's-', 'Color', col_lump, 'LineWidth', lw, ...
         'MarkerFaceColor', col_lump, 'MarkerSize', ms, 'DisplayName', lbl_L);
end
if has_dns
    plot(ax(1), x_bounds_DNS*1000, results_DNS.f_air_seg, '^-', 'Color', col_dns, 'LineWidth', lw_dns, ...
         'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, 'DisplayName', lbl_DNS);
end
plot(ax(1), x_bounds_D*1000, results_D.f_air_seg, 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(1), 'Location','southeast');

if ~isfield(figs_in,'heat_transfer')
    xlabel(ax(2),'Position along HX [mm]'); ylabel(ax(2),'Friction factor f [-]');
    title(ax(2),'Air-side Zanke friction factor through HX');
end
if has_lumped
    plot(ax(2), x_bounds_L*1000, [results_L.f_hx, results_L.f_hx], 's-', 'Color', col_lump, 'LineWidth', lw, ...
         'MarkerFaceColor', col_lump, 'MarkerSize', ms, 'DisplayName', lbl_L);
end
plot(ax(2), x_bounds_D*1000, results_D.f_hx_seg, 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(2), 'Location','southeast');

if ~isfield(figs_in,'heat_transfer')
    xlabel(ax(3),'Position along HX [mm]'); ylabel(ax(3),'Nusselt number [-]');
    title(ax(3),'Nusselt number through HX');
end
if has_lumped
    plot(ax(3), x_bounds_L*1000, [results_L.Nu_air, results_L.Nu_air], 's-', 'Color', col_lump, 'LineWidth', lw, ...
         'MarkerFaceColor', col_lump, 'MarkerSize', ms, 'DisplayName', lbl_L);
end
if has_dns
    plot(ax(3), x_mid_DNS*1000, results_DNS.Nu_air_seg, '^-', 'Color', col_dns, 'LineWidth', lw_dns, ...
         'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, 'DisplayName', lbl_DNS);
end
plot(ax(3), x_mid_D*1000, results_D.Nu_air_seg, 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(3), 'Location','northeast');
figs_out.heat_transfer = fig;

%% ---- Figure: Thermal Performance (h, q, cumulative Q) ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'thermal_perf', 3, 1, 'HX Thermal Performance');

if ~isfield(figs_in,'thermal_perf')
    xlabel(ax(1),'Position along HX [mm]'); ylabel(ax(1),'h [W/m^2K]');
    title(ax(1),'Air-side heat transfer coefficient through HX');
end
if has_lumped
    plot(ax(1), x_bounds_L*1000, [results_L.h_air, results_L.h_air], 's-', 'Color', col_lump, 'LineWidth', lw, ...
         'MarkerFaceColor', col_lump, 'MarkerSize', ms, 'DisplayName', lbl_L);
end
if has_dns
    plot(ax(1), x_mid_DNS*1000, results_DNS.h_air_seg, '^-', 'Color', col_dns, 'LineWidth', lw_dns, ...
         'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, 'DisplayName', lbl_DNS);
end
plot(ax(1), x_mid_D*1000, results_D.h_air_seg, 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(1), 'Location','northeast');

if ~isfield(figs_in,'thermal_perf')
    xlabel(ax(2),'Position along HX [mm]'); ylabel(ax(2),'q [kW/m^2]');
    title(ax(2),'Air-side heat flux through HX');
end
if has_lumped
    plot(ax(2), x_bounds_L*1000, [results_L.h_air, results_L.h_air]*(mean([T_h_i, results_L.T_cool_out]) - results_L.T_mean_c)/1000, ...
         's-', 'Color', col_lump, 'LineWidth', lw, 'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
         'DisplayName', lbl_L);
end
if has_dns
    mean_T_airs_DNS = movmean(results_DNS.T_air_seg, 2); mean_T_airs_DNS = mean_T_airs_DNS(2:end);
    mean_T_cools_DNS = (T_h_i + results_DNS.T_cool_out_arr)/2;
    plot(ax(2), x_mid_DNS*1000, results_DNS.h_air_seg.*(mean_T_cools_DNS - mean_T_airs_DNS)/1000, '^-', ...
         'Color', col_dns, 'LineWidth', lw_dns, 'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, ...
         'DisplayName', lbl_DNS);
end
mean_T_airs_D = movmean(results_D.T_air_seg, 2); mean_T_airs_D = mean_T_airs_D(2:end);
mean_T_cools_D = (T_h_i + results_D.T_cool_out_arr)/2;
plot(ax(2), x_mid_D*1000, results_D.h_air_seg.*(mean_T_cools_D - mean_T_airs_D)/1000, 'o-', ...
     'Color', col_disc, 'LineWidth', lw_disc, 'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(2), 'Location','northeast');

if ~isfield(figs_in,'thermal_perf')
    xlabel(ax(3),'Position along HX [mm]'); ylabel(ax(3),'\Delta Q [kW]');
    title(ax(3),'Cumulative heat rejection through HX');
end
if has_lumped
    plot(ax(3), x_bounds_L*1000, [0, results_L.Q_pred_solution]/1000, 's-', 'Color', col_lump, 'LineWidth', lw, ...
         'MarkerFaceColor', col_lump, 'MarkerSize', ms, 'DisplayName', lbl_L);
end
if has_dns
    plot(ax(3), x_bounds_DNS*1000, [0, cumsum(results_DNS.Q_seg_arr)/1000], '^-', 'Color', col_dns, ...
         'LineWidth', lw_dns, 'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, 'DisplayName', lbl_DNS);
end
plot(ax(3), x_bounds_D*1000, [0, cumsum(results_D.Q_seg_arr)/1000], 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(3), 'Location','northeast');
figs_out.thermal_perf = fig;

%% ---- Figure: Pressure Drop ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'pressure_drop', 2, 1, 'HX Pressure Drop');

if ~isfield(figs_in,'pressure_drop')
    xlabel(ax(1),'Segment midpoint position [mm]'); ylabel(ax(1),'\Delta p per segment [Pa]');
    title(ax(1),'Pressure drop per segment');
end
plot(ax(1), x_mid_D*1000, results_D.dp_seg_arr, 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(1), 'Location','northeast');

if ~isfield(figs_in,'pressure_drop')
    xlabel(ax(2),'Position along HX length [mm]'); ylabel(ax(2),'Cumulative \Delta p [Pa]');
    title(ax(2),'Cumulative pressure drop through HX');
end
if has_lumped
    plot(ax(2), x_bounds_L*1000, [0, results_L.dp_hx], 's-', 'Color', col_lump, 'LineWidth', lw, ...
         'MarkerFaceColor', col_lump, 'MarkerSize', ms, 'DisplayName', lbl_L);
end
if has_dns
    plot(ax(2), x_bounds_DNS*1000, [0, cumsum(results_DNS.dp_seg_arr)], '^-', 'Color', col_dns, ...
         'LineWidth', lw_dns, 'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, 'DisplayName', lbl_DNS);
end
plot(ax(2), x_bounds_D*1000, [0, cumsum(results_D.dp_seg_arr)], 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(2), 'Location','northwest');
figs_out.pressure_drop = fig;

%% ---- Figure: Mass Flow Conservation ---- %%
[fig, ax] = get_or_create_fig(figs_in, 'mass_flow', 1, 1, 'Mass Flow Conservation');

if ~isfield(figs_in,'mass_flow')
    % yline(ax(1), M_dot_3, '--k', 'LineWidth', lw, 'DisplayName', 'M\_dot\_3 (HX inlet)', 'LabelHorizontalAlignment','left');
    % yline(ax(1), M_dot_4, '--b', 'LineWidth', lw, 'DisplayName', 'M\_dot\_4 (HX outlet)', 'LabelHorizontalAlignment','right');
    xlabel(ax(1),'Position along HX [mm]'); ylabel(ax(1),'Air mass flow rate [kg/s]');
    title(ax(1),'Air mass flow rate conservation through HX');
end

if has_dns
    plot(ax(1), x_bounds_DNS*1000, results_DNS.m_dot_seg, '^-', 'Color', col_dns, 'LineWidth', lw_dns, ...
         'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, 'DisplayName', lbl_DNS);
end
plot(ax(1), x_bounds_D*1000, results_D.m_dot_seg, 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(1), 'Location','best');
figs_out.mass_flow = fig;

fprintf("\n=== MASS FLOW CONSERVATION (N=%d) ===\n", N_segments);
fprintf("M_dot_3             = %.6f kg/s\n", M_dot_3);
fprintf("M_dot_4             = %.6f kg/s\n", M_dot_4);
fprintf("Max deviation       = %.2e kg/s\n", max(abs(results_D.m_dot_seg - M_dot_3)));


%% ---- Figure: Boundary layer evolution ---- %%

[fig, ax] = get_or_create_fig(figs_in, 'boundary_layer', 2, 1, 'Boundary Layer Analysis');

if ~isfield(figs_in,'boundary_layer')
    xlabel(ax(1),'Segment midpoint position [mm]'); ylabel(ax(1),'\Boundary layer height [mm]');
    title(ax(1),'Boundary layer height through HX');
end
if has_dns
    plot(ax(1), x_mid_DNS*1000, results_DNS.delta_BL_seg*1000, '^-', 'Color', col_dns, ...
        'LineWidth', lw_dns, 'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, 'DisplayName', lbl_DNS);
end
plot(ax(1), x_mid_D*1000, results_D.delta_BL_seg*1000, 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
    'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(1), 'Location','northeast');

if ~isfield(figs_in,'boundary_layer')
    xlabel(ax(2),'Segment midpoint position [mm]'); ylabel(ax(2),'Hydraulic diameter [mm]');
    title(ax(2),'Hydraulic diameter through HX');
end
if has_dns
    plot(ax(2), x_mid_DNS*1000, results_DNS.d_h_bulk_seg*1000, '^-', 'Color', col_dns, ...
        'LineWidth', lw_dns, 'MarkerFaceColor', col_dns, 'MarkerSize', ms_dns, 'DisplayName', lbl_DNS);
end
plot(ax(2), x_mid_D*1000, results_D.d_h_bulk_seg*1000, 'o-', 'Color', col_disc, 'LineWidth', lw_disc, ...
    'MarkerFaceColor', col_disc, 'MarkerSize', ms_disc, 'DisplayName', lbl);
legend(ax(2), 'Location','northwest');
figs_out.pressure_drop = fig;

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





% function plot_HX(results_D, varargin)
% % Usage:
% %   plot_HX(results_D)                          — discretised only
% %   plot_HX(results_D, results_L)               — + lumped comparison
% %   plot_HX(results_D, results_L, results_DNS)  — + DNS comparison
% %   plot_HX(results_D, [], results_DNS)         — DNS only, no lumped
% 
% has_lumped = ~isempty(varargin) && ~isempty(varargin{1});
% has_dns    = numel(varargin) >= 2 && ~isempty(varargin{2});
% 
% if has_lumped
%     results_L = varargin{1};
% end
% if has_dns
%     results_DNS = varargin{2};
% end
% 
% %% ---- Unpack discretised results ---- %%
% 
% N_segments    = results_D.N_segments;
% L_solution_D  = results_D.L_solution;
% T_h_i         = results_D.T_h_i;
% T_c_i         = results_D.T_c_i;
% T_air_seg     = results_D.T_air_seg;
% P_air_seg     = results_D.P_air_seg;
% v_air_seg     = results_D.v_air_seg;
% Re_air_seg    = results_D.Re_air_seg;
% Pr_air_seg    = results_D.Pr_air_seg;
% f_air_seg     = results_D.f_air_seg;
% Nu_air_seg    = results_D.Nu_air_seg;
% h_air_seg     = results_D.h_air_seg;
% Q_seg_arr     = results_D.Q_seg_arr;
% dp_seg_arr    = results_D.dp_seg_arr;
% T_cool_out_arr = results_D.T_cool_out_arr;
% dp_cool_seg   = results_D.dp_cool_seg;
% v_channel_seg = results_D.v_channel_seg;
% P3_D          = results_D.P3;
% P4_D          = results_D.P4;
% M_dot_3       = results_D.M_dot_3;
% M_dot_4       = results_D.M_dot_4;
% dp_hx_D       = results_D.dp_hx;
% m_dot_seg     = results_D.m_dot_seg;
% f_hx_seg      = results_D.f_hx_seg;
% 
% %% ---- Unpack lumped results ---- %%
% 
% if has_lumped
%     L_solution_L = results_L.L_solution;
%     T3_L         = results_L.T3;
%     T4_L         = results_L.T4;
%     P3_L         = results_L.P3;
%     P4_L         = results_L.P4;
%     dp_hx_L      = results_L.dp_hx;
%     h_air_L      = results_L.h_air;
%     Q_L          = results_L.Q_pred_solution;
%     T_h_o_L      = results_L.T_cool_out;
%     dp_cool_hx_L = results_L.dp_coolant;
%     v_channel_L  = results_L.v_channel_air;
%     Re_air_L     = results_L.Re_air;
%     Pr_air_L     = results_L.Pr_air;
%     f_air_L      = results_L.f_air;
%     Nu_air_L     = results_L.Nu_air;
%     T_mean_c_L   = results_L.T_mean_c;
%     f_hx_L       = results_L.f_hx;
% end
% 
% %% ---- Unpack DNS-based results ---- %%
% 
% if has_dns
%     N_segments_DNS    = results_DNS.N_segments;
%     L_solution_DNS    = results_DNS.L_solution;
%     T_air_seg_DNS     = results_DNS.T_air_seg;
%     P_air_seg_DNS     = results_DNS.P_air_seg;
%     v_air_seg_DNS     = results_DNS.v_air_seg;
%     Re_air_seg_DNS    = results_DNS.Re_air_seg;
%     Pr_air_seg_DNS    = results_DNS.Pr_air_seg;
%     f_air_seg_DNS     = results_DNS.f_air_seg;
%     Nu_air_seg_DNS    = results_DNS.Nu_air_seg;
%     h_air_seg_DNS     = results_DNS.h_air_seg;
%     Q_seg_arr_DNS     = results_DNS.Q_seg_arr;
%     dp_seg_arr_DNS    = results_DNS.dp_seg_arr;
%     T_cool_out_arr_DNS = results_DNS.T_cool_out_arr;
%     dp_cool_seg_DNS   = results_DNS.dp_cool_seg;
%     v_channel_seg_DNS = results_DNS.v_channel_seg;
%     dp_hx_DNS         = results_DNS.dp_hx;
%     m_dot_seg_DNS     = results_DNS.m_dot_seg;
% end
% 
% %% ---- Axis positions ---- %%
% 
% x_bounds_D = linspace(0, L_solution_D, N_segments+1);
% x_mid_D    = 0.5*(x_bounds_D(1:end-1) + x_bounds_D(2:end));
% 
% if has_lumped
%     x_bounds_L = [0, L_solution_L];
% end
% if has_dns
%     x_bounds_DNS = linspace(0, L_solution_DNS, N_segments_DNS+1);
%     x_mid_DNS    = 0.5*(x_bounds_DNS(1:end-1) + x_bounds_DNS(2:end));
% end
% 
% %% ---- Colour scheme ---- %%
% col_lump = [0.122 0.471 0.706];   % blue   — lumped
% col_disc = [0.839 0.153 0.157];   % red    — discretised
% col_dns  = [0.580 0.404 0.741];   % purple — DNS-based
% lw = 1.5;
% ms = 6;
% 
% %% ---- Figure 1: Air temperature and pressure ---- %%
% figure('Name','HX Air State Variables','NumberTitle','off');
% 
% subplot(2,1,1)
% hold on
% plot(x_bounds_D*1000, T_air_seg, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
%      'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
%      'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [T3_L, T4_L], 's-', 'Color', col_lump, 'LineWidth', lw, ...
%          'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%          'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
%     plot(x_bounds_L*1000, [T_mean_c_L, T_mean_c_L], 's--', 'Color', col_lump, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_lump, 'MarkerSize', ms, 'DisplayName', 'Lumped T_{mean,c}')
% end
% if has_dns
%     plot(x_bounds_DNS*1000, T_air_seg_DNS, '^-', 'Color', col_dns, 'LineWidth', lw, ...
%          'MarkerFaceColor', col_dns, 'MarkerSize', ms, ...
%          'DisplayName', sprintf('DNS-based (N=%d, L=%.3fm)', N_segments_DNS, L_solution_DNS))
% end
% yline(T_h_i, '--r', 'DisplayName', 'T_{h,i}', 'LabelHorizontalAlignment','left')
% yline(T_c_i, '--b', 'DisplayName', 'T_{c,i}', 'LabelHorizontalAlignment','left')
% xlabel('Position along HX [mm]')
% ylabel('Air temperature [K]')
% title('Air temperature distribution through HX')
% legend('Location','southeast')
% grid on
% 
% subplot(2,1,2)
% hold on
% plot(x_bounds_D*1000, P_air_seg/1000, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
%      'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
%      'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [P3_L, P4_L]/1000, 's-', 'Color', col_lump, 'LineWidth', lw, ...
%          'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%          'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
% end
% if has_dns
%     plot(x_bounds_DNS*1000, P_air_seg_DNS/1000, '^-', 'Color', col_dns, 'LineWidth', lw, ...
%          'MarkerFaceColor', col_dns, 'MarkerSize', ms, ...
%          'DisplayName', sprintf('DNS-based (N=%d, L=%.3fm)', N_segments_DNS, L_solution_DNS))
% end
% xlabel('Position along HX [mm]')
% ylabel('Air static pressure [kPa]')
% title('Air pressure distribution through HX')
% legend('Location','northeast')
% grid on
% 
% %% ---- Figure: Coolant temperature and pressure drop ---- %%
% figure('Name','HX Coolant Variables','NumberTitle','off');
% 
% subplot(2,1,1)
% hold on
% plot(x_mid_D*1000, T_cool_out_arr, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
%     'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
%     'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [T_h_o_L, T_h_o_L], 's-', 'Color', col_lump, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
% end
% if has_dns
%     plot(x_mid_DNS*1000, T_cool_out_arr_DNS, '^-', 'Color', col_dns, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_dns, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('DNS-based (N=%d, L=%.3fm)', N_segments_DNS, L_solution_DNS))
% end
% yline(T_h_i, '--r', 'DisplayName', 'T_{h,i}', 'LabelHorizontalAlignment','left')
% yline(T_c_i, '--b', 'DisplayName', 'T_{c,i}', 'LabelHorizontalAlignment','left')
% xlabel('Position along HX [mm]')
% ylabel('Coolant temperature [K]')
% title('Coolant temperature distribution through HX')
% legend('Location','southeast')
% grid on
% 
% subplot(2,1,2)
% hold on
% plot(x_mid_D*1000, dp_cool_seg/1e5, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
%     'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
%     'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [dp_cool_hx_L, dp_cool_hx_L]/1e5, 's-', 'Color', col_lump, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
% end
% if has_dns
%     plot(x_mid_DNS*1000, dp_cool_seg_DNS/1e5, '^-', 'Color', col_dns, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_dns, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('DNS-based (N=%d, L=%.3fm)', N_segments_DNS, L_solution_DNS))
% end
% xlabel('Position along HX [mm]')
% ylabel('Coolant pressure drop [bar]')
% title('Coolant pressure drop distribution through HX')
% legend('Location','southeast')
% grid on
% 
% %% ---- Figure 2: Air velocity and Reynolds number ---- %%
% figure('Name','HX Air Flow Variables','NumberTitle','off');
% 
% subplot(3,1,1)
% hold on
% plot(x_bounds_D*1000, v_channel_seg, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
%     'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
%     'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [v_channel_L, v_channel_L], 's-', 'Color', col_lump, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
% end
% if has_dns
%     plot(x_bounds_DNS*1000, v_channel_seg_DNS, '^-', 'Color', col_dns, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_dns, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('DNS-based (N=%d, L=%.3fm)', N_segments_DNS, L_solution_DNS))
% end
% xlabel('Position along HX [mm]')
% ylabel('Air channel velocity [m/s]')
% title('Channel air velocity distribution through HX')
% legend('location', 'northeast');
% grid on
% 
% subplot(3,1,2)
% hold on
% plot(x_bounds_D*1000, Re_air_seg, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
%     'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
%     'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [Re_air_L, Re_air_L], 's-', 'Color', col_lump, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
% end
% if has_dns
%     plot(x_bounds_DNS*1000, Re_air_seg_DNS, '^-', 'Color', col_dns, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_dns, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('DNS-based (N=%d, L=%.3fm)', N_segments_DNS, L_solution_DNS))
% end
% xlabel('Position along HX [mm]')
% ylabel('Reynolds number [-]')
% title('Air Reynolds number through HX')
% yline(2300, '--k', 'DisplayName', 'Re = 2300', 'LabelHorizontalAlignment','left')
% yline(3000, '--r', 'DisplayName', 'Re = 3000', 'LabelHorizontalAlignment','left')
% legend('location', 'northeast');
% grid on
% 
% subplot(3,1,3)
% hold on
% plot(x_bounds_D*1000, Pr_air_seg, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
%     'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
%     'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [Pr_air_L, Pr_air_L], 's-', 'Color', col_lump, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
% end
% if has_dns
%     plot(x_bounds_DNS*1000, Pr_air_seg_DNS, '^-', 'Color', col_dns, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_dns, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('DNS-based (N=%d, L=%.3fm)', N_segments_DNS, L_solution_DNS))
% end
% xlabel('Position along HX [mm]')
% ylabel('Prandtl number [-]')
% title('Air Prandtl number through HX')
% legend('location', 'northeast');
% grid on
% 
% %% ---- Figure 3: Friction factor and Nusselt number ---- %%
% figure('Name','HX Heat Transfer Variables','NumberTitle','off');
% 
% subplot(3,1,1)
% hold on
% plot(x_bounds_D*1000, f_air_seg, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
%     'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
%     'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [f_air_L, f_air_L], 's-', 'Color', col_lump, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
% end
% if has_dns
%     plot(x_bounds_DNS*1000, f_air_seg_DNS, '^-', 'Color', col_dns, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_dns, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('DNS-based (N=%d, L=%.3fm)', N_segments_DNS, L_solution_DNS))
% end
% xlabel('Position along HX [mm]')
% ylabel('Friction factor f [-]')
% title('Air-side Konakov friction factor through HX')
% legend('location', 'southeast');
% grid on
% 
% subplot(3,1,2)
% hold on
% plot(x_bounds_D*1000, f_hx_seg, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
%     'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
%     'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [f_hx_L, f_hx_L], 's-', 'Color', col_lump, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
% end
% xlabel('Position along HX [mm]')
% ylabel('Friction factor f [-]')
% title('Air-side Zanke friction factor through HX')
% legend('location', 'southeast');
% grid on
% 
% subplot(3,1,3)
% hold on
% plot(x_mid_D*1000, Nu_air_seg, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
%     'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
%     'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [Nu_air_L, Nu_air_L], 's-', 'Color', col_lump, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
% end
% if has_dns
%     plot(x_mid_DNS*1000, Nu_air_seg_DNS, '^-', 'Color', col_dns, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_dns, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('DNS-based (N=%d, L=%.3fm)', N_segments_DNS, L_solution_DNS))
% end
% xlabel('Position along HX [mm]')
% ylabel('Nusselt number [-]')
% title('Nusselt number through HX')
% legend('location', 'northeast');
% grid on
% 
% %% ---- Figure 4: Heat transfer coefficient and heat rejected ---- %%
% figure('Name','HX Thermal Performance','NumberTitle','off');
% 
% subplot(3,1,1)
% hold on
% plot(x_mid_D*1000, h_air_seg, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
%     'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
%     'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [h_air_L, h_air_L], 's-', 'Color', col_lump, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
% end
% if has_dns
%     plot(x_mid_DNS*1000, h_air_seg_DNS, '^-', 'Color', col_dns, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_dns, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('DNS-based (N=%d, L=%.3fm)', N_segments_DNS, L_solution_DNS))
% end
% xlabel('Position along HX [mm]')
% ylabel('h [W/m^2K]')
% title('Air-side heat transfer coefficient through HX')
% legend('location', 'northeast');
% grid on
% 
% subplot(3,1,2)
% hold on
% mean_T_airs = movmean(T_air_seg, 2);
% mean_T_airs = mean_T_airs(2:end);
% mean_T_cools = (T_h_i+T_cool_out_arr)/2;
% 
% plot(x_mid_D*1000, h_air_seg.*(mean_T_cools - mean_T_airs)/1000, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
%     'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
%     'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [h_air_L, h_air_L]*(mean([T_h_i, T_h_o_L]) - T_mean_c_L)/1000, 's-', 'Color', col_lump, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
% end
% if has_dns
%     mean_T_airs_DNS = movmean(T_air_seg_DNS, 2);
%     mean_T_airs_DNS = mean_T_airs_DNS(2:end);
%     mean_T_cools_DNS = (T_h_i+T_cool_out_arr_DNS)/2;
%     plot(x_mid_DNS*1000, h_air_seg_DNS.*(mean_T_cools_DNS - mean_T_airs_DNS)/1000, '^-', 'Color', col_dns, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_dns, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('DNS-based (N=%d, L=%.3fm)', N_segments_DNS, L_solution_DNS))
% end
% xlabel('Position along HX [mm]')
% ylabel('q [kW/m^2]')
% title('Air-side heat flux through HX')
% legend('location', 'northeast');
% grid on
% 
% subplot(3,1,3)
% hold on
% plot(x_bounds_D*1000, [0, cumsum(Q_seg_arr)/1000], 'o-', 'Color', col_disc, 'LineWidth', lw, ...
%     'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
%     'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [0, Q_L]/1000, 's-', 'Color', col_lump, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
% end
% if has_dns
%     plot(x_bounds_DNS*1000, [0, cumsum(Q_seg_arr_DNS)/1000], '^-', 'Color', col_dns, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_dns, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('DNS-based (N=%d, L=%.3fm)', N_segments_DNS, L_solution_DNS))
% end
% xlabel('Position along HX [mm]')
% ylabel('\Delta Q [kW]')
% title('Cumulative heat rejection through HX')
% yline(sum(Q_seg_arr)/N_segments/1000, '--k', 'DisplayName', 'Mean Q/segment', ...
%     'LabelHorizontalAlignment','left')
% legend('location', 'northeast');
% grid on
% 
% %% ---- Figure 5: Pressure drop ---- %%
% figure('Name','HX Pressure Drop','NumberTitle','off');
% 
% subplot(2,1,1)
% hold on
% bar(x_mid_D*1000, dp_seg_arr, 'FaceColor', col_disc, 'EdgeColor', 'none')
% xlabel('Segment midpoint position [mm]')
% ylabel('\Delta p per segment [Pa]')
% title(sprintf('Pressure drop per segment (N=%d)', N_segments))
% grid on
% 
% subplot(2,1,2)
% hold on
% plot(x_bounds_D*1000, [0, cumsum(dp_seg_arr)], 'o-', 'Color', col_disc, ...
%      'LineWidth', lw, 'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
%      'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [0, dp_hx_L], 's-', 'Color', col_lump, ...
%          'LineWidth', lw, 'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%          'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
% end
% if has_dns
%     plot(x_bounds_DNS*1000, [0, cumsum(dp_seg_arr_DNS)], '^-', 'Color', col_dns, ...
%          'LineWidth', lw, 'MarkerFaceColor', col_dns, 'MarkerSize', ms, ...
%          'DisplayName', sprintf('DNS-based (N=%d, L=%.3fm)', N_segments_DNS, L_solution_DNS))
% end
% xlabel('Position along HX length [mm]')
% ylabel('Cumulative \Delta p [Pa]')
% title('Cumulative pressure drop through HX')
% legend('Location','northwest')
% grid on
% 
% %% ---- Figure 6: Mass flow conservation ---- %%
% 
% figure('Name','Mass Flow Conservation','NumberTitle','off');
% hold on
% plot(x_bounds_D*1000, m_dot_seg, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
%      'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
%      'DisplayName', sprintf('Discretised (N=%d)', N_segments))
% if has_dns
%     plot(x_bounds_DNS*1000, m_dot_seg_DNS, '^-', 'Color', col_dns, 'LineWidth', lw, ...
%          'MarkerFaceColor', col_dns, 'MarkerSize', ms, ...
%          'DisplayName', sprintf('DNS-based (N=%d)', N_segments_DNS))
% end
% yline(M_dot_3, '--k', 'LineWidth', lw, 'DisplayName', 'M\_dot\_3 (HX inlet)', 'LabelHorizontalAlignment','left')
% yline(M_dot_4, '--b', 'LineWidth', lw, 'DisplayName', 'M\_dot\_4 (HX outlet)', 'LabelHorizontalAlignment','right')
% xlabel('Position along HX [mm]')
% ylabel('Air mass flow rate [kg/s]')
% title('Air mass flow rate conservation through HX')
% legend('Location','best')
% grid on
% 
% fprintf("\n=== MASS FLOW CONSERVATION ===\n");
% fprintf("M_dot_3             = %.6f kg/s\n", M_dot_3);
% fprintf("M_dot_4             = %.6f kg/s\n", M_dot_4);
% fprintf("Max deviation       = %.2e kg/s\n", max(abs(m_dot_seg - M_dot_3)));
% 
% end