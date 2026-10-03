function comparison = compare_cfd_vs_2d(results_2D, results_cfd, varargin)
% compare_cfd_vs_2d  Plot CFD unit-cell results (from process_wall_CFD_data)
% against a single coolant-side segment of the 2D discretised model - by
% default the FIRST one (l=1) - each at its own native x positions, with
% no interpolation. Compares air bulk temperature, h, Nu, static
% pressure, velocity, and mass flow rate.
%
% Rationale for l=1 (default): the CFD unit cell models coolant crossing
% the air channel once, entering fresh at T_h_i - exactly the model's
% l=1 boundary condition (T_cool_seg(:,1) = T_h_i for every air segment,
% from the inlet-manifold assumption). Segments l>1 represent coolant
% that has already warmed up crossing earlier width-positions, which a
% single-crossing CFD unit cell does not capture - so those are NOT a
% fair comparison target.
%
% Usage:
%   comparison = compare_cfd_vs_2d(results_2D, results_cfd);
%   comparison = compare_cfd_vs_2d(results_2D, results_cfd, 'l_index', 2);
%
% results_2D  : results struct from the 2D discretised model, carrying
%               T_air_seg, h_air_seg, Nu_air_seg, P_air_seg,
%               v_channel_seg, m_dot_seg (N_segments(+1) x N_seg_cool),
%               L_solution, N_segments.
% results_cfd : struct returned by process_wall_CFD_data (z, Tbulk, h, Nu,
%               dP_bulk, Vmag, mdot, is_entrance).
%
% Entrance (z=0) station: only h/Nu have a genuine thermal-entrance
% singularity there (h formally diverges for a step-change wall
% condition), so only their panels exclude it. Temperature, pressure,
% velocity and mass flow are all well-behaved at z=0 and include it - for
% temperature/velocity/mass flow it's the actual inlet condition; for
% pressure it's the anchor point where the two traces coincide by
% construction (see below).
%
% Pressure: results_cfd.dP_bulk is CFD's static pressure already shifted
% so its z=0 (inlet) station reads zero (see process_wall_CFD_data) -
% adding the model's own inlet static pressure, P_air_seg(1,l), gives an
% absolute CFD pressure trace directly comparable to the model's
% P_air_seg(:,l).
%
% No interpolation: the CFD z-grid is fixed by the mesh and the model's
% node/segment positions depend on the SOLVED L_solution, so the two
% series are plotted at their own native x stations, on shared axes, and
% left to visually align (or not) rather than being resampled onto a
% common grid.

p = inputParser;
addParameter(p, 'l_index', 1);
parse(p, varargin{:});
l = p.Results.l_index;

if l > size(results_2D.T_air_seg, 2)
    error('compare_cfd_vs_2d:lOutOfRange', ...
        'l_index = %d requested but results_2D only has %d coolant-side segments.', ...
        l, size(results_2D.T_air_seg, 2));
end

%% ---- Model x-positions [mm] ---- %%
x_bounds = linspace(0, results_2D.L_solution, results_2D.N_segments+1) * 1000;  % node positions
x_mid    = 0.5*(x_bounds(1:end-1) + x_bounds(2:end));                            % segment midpoints

T_model  = results_2D.T_air_seg(:, l)';
h_model  = results_2D.h_air_seg(:, l)';
Nu_model = results_2D.Nu_air_seg(:, l)';
P_model  = results_2D.P_air_seg(:, l)';
v_model  = results_2D.v_channel_seg(:, l)';
mdot_model = results_2D.m_dot_seg./results_2D.N_fin_air*2./results_2D.N_air_pass;

%% ---- CFD, all stations: T, P, V and mdot are non-singular at z=0 (inlet) ---- %%
x_cfd_all  = -results_cfd.z * 1000;        % flip sign (flow runs in -z), m -> mm
Tb_cfd_all = results_cfd.Tbulk;
V_cfd_all  = results_cfd.Vmag;
mdot_cfd_all = results_cfd.mdot;

%% ---- CFD, entrance excluded: h/Nu only (thermal-entrance singularity) ---- %%
valid  = ~results_cfd.is_entrance;
x_cfd  = -results_cfd.z(valid) * 1000;
h_cfd  = results_cfd.h(valid);
Nu_cfd = results_cfd.Nu(valid);

%% ---- CFD static pressure: shifted trace + model's own inlet pressure ---- %%
P_2D_inlet = results_2D.P_air_seg(1, l);
P_cfd = P_2D_inlet + results_cfd.dP_bulk;

col_model = [0.839 0.153 0.157];  % red  - matches plot_HX's 'discretised'/'2D' family
col_cfd   = [0.3 0.3 0.3];        % gray - matches plot_HX's 'cfd' style
lw = 1.2; ms = 5;

%% ---- Figure: air bulk temperature ---- %%
figure('Name', sprintf('Air Bulk Temperature: Model (l=%d) vs CFD', l), 'NumberTitle', 'off');
hold on; grid on;
plot(x_bounds, T_model, 'o-', 'Color', col_model, 'LineWidth', lw, ...
    'MarkerFaceColor', col_model, 'MarkerSize', ms, 'DisplayName', sprintf('2D model (l=%d)', l));
plot(x_cfd_all, Tb_cfd_all, 'x-', 'Color', col_cfd, 'LineWidth', lw, ...
    'MarkerSize', ms+2, 'DisplayName', 'CFD');
xlabel('Position along HX [mm]'); ylabel('Air bulk temperature [K]');
title('Air bulk temperature: model vs CFD');
legend('Location', 'best');

%% ---- Figure: static pressure ---- %%
figure('Name', sprintf('Static Pressure: Model (l=%d) vs CFD', l), 'NumberTitle', 'off');
hold on; grid on;
plot(x_bounds, P_model, 'o-', 'Color', col_model, 'LineWidth', lw, ...
    'MarkerFaceColor', col_model, 'MarkerSize', ms, 'DisplayName', sprintf('2D model (l=%d)', l));
plot(x_cfd_all, P_cfd, 'x-', 'Color', col_cfd, 'LineWidth', lw, ...
    'MarkerSize', ms+2, 'DisplayName', 'CFD (inlet-referenced)');
xlabel('Position along HX [mm]'); ylabel('Static pressure [Pa]');
title('Static pressure: model vs CFD');
legend('Location', 'best');

%% ---- Figure: heat transfer coefficient ---- %%
figure('Name', sprintf('Heat Transfer Coefficient: Model (l=%d) vs CFD', l), 'NumberTitle', 'off');
hold on; grid on;
plot(x_mid, h_model, 'o-', 'Color', col_model, 'LineWidth', lw, ...
    'MarkerFaceColor', col_model, 'MarkerSize', ms, 'DisplayName', sprintf('2D model (l=%d)', l));
plot(x_cfd, h_cfd, 'x-', 'Color', col_cfd, 'LineWidth', lw, ...
    'MarkerSize', ms+2, 'DisplayName', 'CFD');
xlabel('Position along HX [mm]'); ylabel('h [W/m^2K]');
title('Heat transfer coefficient: model vs CFD');
legend('Location', 'best');

%% ---- Figure: velocity ---- %%
figure('Name', sprintf('Velocity: Model (l=%d) vs CFD', l), 'NumberTitle', 'off');
hold on; grid on;
plot(x_bounds, v_model, 'o-', 'Color', col_model, 'LineWidth', lw, ...
    'MarkerFaceColor', col_model, 'MarkerSize', ms, 'DisplayName', sprintf('2D model (l=%d)', l));
plot(x_cfd_all, V_cfd_all, 'x-', 'Color', col_cfd, 'LineWidth', lw, ...
    'MarkerSize', ms+2, 'DisplayName', 'CFD');
xlabel('Position along HX [mm]'); ylabel('Velocity [m/s]');
title('Channel velocity: model vs CFD');
legend('Location', 'best');

%% ---- Figure: mass flow rate ---- %%
figure('Name', sprintf('Mass Flow Rate: Model (l=%d) vs CFD', l), 'NumberTitle', 'off');
hold on; grid on;
plot(x_bounds, mdot_model, 'o-', 'Color', col_model, 'LineWidth', lw, ...
    'MarkerFaceColor', col_model, 'MarkerSize', ms, 'DisplayName', sprintf('2D model (l=%d)', l));
plot(x_cfd_all, mdot_cfd_all, 'x-', 'Color', col_cfd, 'LineWidth', lw, ...
    'MarkerSize', ms+2, 'DisplayName', 'CFD');
xlabel('Position along HX [mm]'); ylabel('Mass flow rate [kg/s]');
title('Mass flow rate: model vs CFD');
legend('Location', 'best');

%% ---- Figure: Nusselt number ---- %%
figure('Name', sprintf('Nusselt Number: Model (l=%d) vs CFD', l), 'NumberTitle', 'off');
hold on; grid on;
plot(x_mid, Nu_model, 'o-', 'Color', col_model, 'LineWidth', lw, ...
    'MarkerFaceColor', col_model, 'MarkerSize', ms, 'DisplayName', sprintf('2D model (l=%d)', l));
plot(x_cfd, Nu_cfd, 'x-', 'Color', col_cfd, 'LineWidth', lw, ...
    'MarkerSize', ms+2, 'DisplayName', 'CFD');
xlabel('Position along HX [mm]'); ylabel('Nu [-]');
title('Nusselt number: model vs CFD');
legend('Location', 'best');

%% ---- Package for programmatic use ---- %%
comparison.l_index    = l;
comparison.x_bounds   = x_bounds;
comparison.x_mid      = x_mid;
comparison.T_model    = T_model;
comparison.h_model    = h_model;
comparison.Nu_model   = Nu_model;
comparison.P_model    = P_model;
comparison.v_model    = v_model;
comparison.mdot_model = mdot_model;
comparison.x_cfd_all  = x_cfd_all;
comparison.T_cfd      = Tb_cfd_all;
comparison.P_cfd      = P_cfd;
comparison.V_cfd      = V_cfd_all;
comparison.mdot_cfd   = mdot_cfd_all;
comparison.x_cfd      = x_cfd;
comparison.h_cfd      = h_cfd;
comparison.Nu_cfd     = Nu_cfd;

end