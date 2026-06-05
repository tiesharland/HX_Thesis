function plot_HX(N_segments, L_solution_D, T_h_i, T_c_i, ...
    T_air_seg, P_air_seg, v_air_seg, Re_air_seg, ...
    f_air_seg, Nu_air_seg, h_air_seg, Q_seg_arr, dp_seg_arr, ...
    T_cool_out_arr, dp_cool_seg, ...
    P3_D, P4_D, M_dot_3, M_dot_4, dp_hx_D, m_dot_seg, ...
    varargin)

has_lumped = ~isempty(varargin);

if has_lumped
    L_solution_L = varargin{1};
    T4_L         = varargin{2};
    P3_L         = varargin{3};
    P4_L         = varargin{4};
    T3_L         = varargin{5};
    dp_hx_L      = varargin{6};
end

%% ---- Axis positions ---- %%

x_bounds_D = linspace(0, L_solution_D, N_segments+1);
x_mid_D    = 0.5*(x_bounds_D(1:end-1) + x_bounds_D(2:end));

if has_lumped
    x_bounds_L = [0, L_solution_L];
end

%% ---- Colour scheme ---- %%
col_air  = [0.122 0.471 0.706];
col_disc = [0.839 0.153 0.157];
col_dp   = [0.173 0.627 0.173];
col_lump = [0.5   0.5   0.5  ];
lw = 1.5;
ms = 6;

%% ---- Figure 1: Air temperature and pressure ---- %%
figure('Name','HX Air State Variables','NumberTitle','off');

subplot(2,1,1)
hold on
plot(x_bounds_D*1000, T_air_seg, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
     'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
if has_lumped
    plot(x_bounds_L*1000, [T3_L, T4_L], 's--', 'Color', col_lump, 'LineWidth', lw, ...
         'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
end
yline(T_h_i, '--r', 'DisplayName', 'T_{h,i}', 'LabelHorizontalAlignment','left')
yline(T_c_i, '--b', 'DisplayName', 'T_{c,i}', 'LabelHorizontalAlignment','left')
xlabel('Position along HX length [mm]')
ylabel('Air temperature [K]')
title('Air temperature distribution through HX')
legend('Location','northwest')
grid on

subplot(2,1,2)
hold on
plot(x_bounds_D*1000, P_air_seg/1000, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
     'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
if has_lumped
    plot(x_bounds_L*1000, [P3_L, P4_L]/1000, 's--', 'Color', col_lump, 'LineWidth', lw, ...
         'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
end
xlabel('Position along HX length [mm]')
ylabel('Air static pressure [kPa]')
title('Air pressure distribution through HX')
legend('Location','northeast')
grid on

%% ---- Figure: Coolant temperature and pressure drop ---- %%
figure('Name','HX Coolant Variables','NumberTitle','off');

subplot(2,1,1)
hold on
plot(x_mid_D*1000, T_cool_out_arr, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
    'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
    'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [T3_L, T4_L], 's--', 'Color', col_lump, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
% end
yline(T_h_i, '--r', 'DisplayName', 'T_{h,i}', 'LabelHorizontalAlignment','left')
yline(T_c_i, '--b', 'DisplayName', 'T_{c,i}', 'LabelHorizontalAlignment','left')
xlabel('Position along HX length [mm]')
ylabel('Coolant temperature [K]')
title('Coolant temperature distribution through HX')
legend('Location','northwest')
grid on

subplot(2,1,2)
hold on
plot(x_mid_D*1000, dp_cool_seg, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
    'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
    'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
% if has_lumped
%     plot(x_bounds_L*1000, [P3_L, P4_L]/1000, 's--', 'Color', col_lump, 'LineWidth', lw, ...
%         'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
%         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
% end
xlabel('Position along HX length [mm]')
ylabel('Coolant pressure drop [Pa]')
title('Coolant pressure drop distribution through HX')
legend('Location','northeast')
grid on

%% ---- Figure 2: Air velocity and Reynolds number ---- %%
figure('Name','HX Air Flow Variables','NumberTitle','off');

subplot(2,1,1)
plot(x_bounds_D*1000, v_air_seg, 'o-', 'Color', col_air, 'LineWidth', lw, ...
     'MarkerFaceColor', col_air, 'MarkerSize', ms)
xlabel('Position along HX length [mm]')
ylabel('Air channel velocity [m/s]')
title(sprintf('Air velocity distribution (N=%d)', N_segments))
grid on

subplot(2,1,2)
plot(x_bounds_D*1000, Re_air_seg, 's-', 'Color', col_air, 'LineWidth', lw, ...
     'MarkerFaceColor', col_air, 'MarkerSize', ms)
xlabel('Position along HX length [mm]')
ylabel('Reynolds number [-]')
title(sprintf('Air Reynolds number (N=%d)', N_segments))
yline(2300, '--k', 'Re = 2300', 'LabelHorizontalAlignment','left')
yline(3000, '--r', 'Re = 3000', 'LabelHorizontalAlignment','left')
grid on

%% ---- Figure 3: Friction factor and Nusselt number ---- %%
figure('Name','HX Heat Transfer Variables','NumberTitle','off');

subplot(2,1,1)
plot(x_bounds_D*1000, f_air_seg, 'o-', 'Color', col_dp, 'LineWidth', lw, ...
     'MarkerFaceColor', col_dp, 'MarkerSize', ms)
xlabel('Position along HX length [mm]')
ylabel('Friction factor f [-]')
title(sprintf('Friction factor per node (N=%d)', N_segments))
grid on

subplot(2,1,2)
plot(x_mid_D*1000, Nu_air_seg, 's-', 'Color', col_air, 'LineWidth', lw, ...
     'MarkerFaceColor', col_air, 'MarkerSize', ms)
xlabel('Segment midpoint position [mm]')
ylabel('Nusselt number [-]')
title(sprintf('Nusselt number per segment (N=%d)', N_segments))
grid on

%% ---- Figure 4: Heat transfer coefficient and heat rejected ---- %%
figure('Name','HX Thermal Performance','NumberTitle','off');

subplot(2,1,1)
plot(x_mid_D*1000, h_air_seg, 'o-', 'Color', col_air, 'LineWidth', lw, ...
     'MarkerFaceColor', col_air, 'MarkerSize', ms)
xlabel('Segment midpoint position [mm]')
ylabel('Heat transfer coefficient [W/m^2K]')
title(sprintf('Air-side h per segment (N=%d)', N_segments))
grid on

subplot(2,1,2)
plot(x_mid_D*1000, Q_seg_arr/1000, 'o-', 'Color', col_disc, 'LineWidth', lw, ...
     'MarkerFaceColor', col_disc, 'MarkerSize', ms)
xlabel('Segment midpoint position [mm]')
ylabel('Heat rejected per segment [kW]')
title(sprintf('Heat rejection per segment (N=%d)', N_segments))
yline(sum(Q_seg_arr)/N_segments/1000, '--k', 'Mean Q/segment', ...
      'LabelHorizontalAlignment','left')
grid on

%% ---- Figure 5: Pressure drop ---- %%
figure('Name','HX Pressure Drop','NumberTitle','off');

subplot(2,1,1)
bar(x_mid_D*1000, dp_seg_arr, 'FaceColor', col_disc, 'EdgeColor', 'none')
xlabel('Segment midpoint position [mm]')
ylabel('\Delta p per segment [Pa]')
title(sprintf('Pressure drop per segment (N=%d)', N_segments))
grid on

subplot(2,1,2)
hold on
plot(x_bounds_D*1000, [0, cumsum(dp_seg_arr)], 'o-', 'Color', col_disc, ...
     'LineWidth', lw, 'MarkerFaceColor', col_disc, 'MarkerSize', ms, ...
     'DisplayName', sprintf('Discretised (N=%d, L=%.3fm)', N_segments, L_solution_D))
if has_lumped
    plot(x_bounds_L*1000, [0, dp_hx_L], 's--', 'Color', col_lump, ...
         'LineWidth', lw, 'MarkerFaceColor', col_lump, 'MarkerSize', ms, ...
         'DisplayName', sprintf('Lumped (L=%.3fm)', L_solution_L))
    yline(dp_hx_L, '--', 'Color', col_lump, 'DisplayName', 'Total dp lumped', ...
          'LabelHorizontalAlignment','right')
end
yline(dp_hx_D, '--', 'Color', col_disc, 'DisplayName', 'Total dp disc', ...
      'LabelHorizontalAlignment','left')
xlabel('Position along HX length [mm]')
ylabel('Cumulative \Delta p [Pa]')
title('Cumulative pressure drop through HX')
legend('Location','northwest')
grid on

%% ---- Figure 6: Mass flow conservation ---- %%

figure('Name','Mass Flow Conservation','NumberTitle','off');
plot(x_bounds_D*1000, m_dot_seg, 'o-', 'Color', col_air, 'LineWidth', lw, ...
     'MarkerFaceColor', col_air, 'MarkerSize', ms)
hold on
yline(M_dot_3, '--k', 'DisplayName', 'M\_dot\_3 (HX inlet)', 'LabelHorizontalAlignment','left')
yline(M_dot_4, '--r', 'DisplayName', 'M\_dot\_4 (HX outlet)', 'LabelHorizontalAlignment','left')
xlabel('Position along HX length [mm]')
ylabel('Air mass flow rate [kg/s]')
title('Air mass flow rate conservation through HX segments')
grid on

fprintf("\n=== MASS FLOW CONSERVATION ===\n");
fprintf("M_dot_3             = %.6f kg/s\n", M_dot_3);
fprintf("M_dot_4             = %.6f kg/s\n", M_dot_4);
fprintf("Max deviation       = %.2e kg/s\n", max(abs(m_dot_seg - M_dot_3)));

end