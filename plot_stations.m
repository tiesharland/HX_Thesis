function plot_stations(P_inf, P_inf_tot, T_inf, M_inf, V_inf, ...
    P1_D, P2_D, P3_D, P4_D, P5_D, P6_D, ...
    T1_D, T2_D, T3_D, T4_D, T5_D, T6_D, ...
    P1_0_D, P2_0_D, P3_0_D, P4_0_D, P5_0_D, P6_0_D, ...
    v1_D, v2_D, v3_D, v4_D, v5_D, v6_D, ...
    M1_D, M2_D, M3_D, M4_D, M5_D, M6_D, ...
    N_segments, T_air_seg, P_air_seg, v_air_seg, v_channel_seg, varargin)

has_lumped = ~isempty(varargin);

if has_lumped
    % % TMS state arrays — lumped
    % P_L   = [P_inf,     P1_L, P2_L, P3_L, P4_L, P5_L, P6_L];
    % P0_L  = [P_inf_tot, P1_0_L, P2_0_L, P3_0_L, P4_0_L, P5_0_L, P6_0_L];
    % T_L   = [T_inf,     T1_L, T2_L, T3_L, T4_L, T5_L, T6_L];
    % v_L   = [V_inf,     v1_L, v2_L, v3_L, v4_L, v5_L, v6_L];
    % M_L   = [M_inf,     M1_L, M2_L, M3_L, M4_L, M5_L, M6_L];

    P_L = [P_inf, varargin{1:6}];
    T_L = [T_inf, varargin{7:12}];
    P0_L = [P_inf_tot, varargin{13:18}];
    v_L = [V_inf, varargin{19:24}];
    M_L = [M_inf, varargin{25:30}];
end

%% ---- Colour scheme ---- %%
col_L  = [0.122 0.471 0.706];   % blue  — lumped
col_D  = [0.839 0.153 0.157];   % red   — discretised
col_HX = [0.173 0.627 0.173];   % green — HX interior (discretised only)
lw     = 1.5;
ms     = 7;

%% ---- Station x-positions ---- %%
% Stations 0 (freestream) through 6
x_stations = 0:6;
labels      = {'0 (∞)','1','2','3','4','5','6'};

% TMS state arrays — discretised
P_D   = [P_inf,     P1_D, P2_D, P3_D, P4_D, P5_D, P6_D];
P0_D  = [P_inf_tot, P1_0_D, P2_0_D, P3_0_D, P4_0_D, P5_0_D, P6_0_D];
T_D   = [T_inf,     T1_D, T2_D, T3_D, T4_D, T5_D, T6_D];
v_D   = [V_inf,     v1_D, v2_D, v3_D, v4_D, v5_D, v6_D];
M_D   = [M_inf,     M1_D, M2_D, M3_D, M4_D, M5_D, M6_D];

% HX interior x-positions for discretised model
% Map node positions (x_bounds) to the station axis between 3 and 4
x_HX = linspace(3, 4, N_segments+1);

%% ---- Figure 1: Static pressure through TMS ---- %%
figure('Name','Static Pressure — TMS Comparison','NumberTitle','off');
hold on

if has_lumped
    plot(x_stations, P_L/1e3, 'o-', 'Color', col_L, 'LineWidth', lw, ...
         'MarkerFaceColor', col_L, 'MarkerSize', ms, 'DisplayName', 'Lumped')
end
plot(x_stations, P_D/1e3, 's--', 'Color', col_D, 'LineWidth', lw, ...
 'MarkerFaceColor', col_D, 'MarkerSize', ms, 'DisplayName', 'Discretised')

% HX interior for discretised
plot(x_HX, P_air_seg/1e3, '-', 'Color', col_HX, 'LineWidth', lw, ...
     'DisplayName', 'Disc. HX interior')

set(gca, 'XTick', x_stations, 'XTickLabel', labels)
xlabel('TMS Station')
ylabel('Static pressure [kPa]')
title('Static pressure through TMS')
legend('Location', 'best')
grid on

%% ---- Figure 2: Total pressure through TMS ---- %%
figure('Name','Total Pressure — TMS Comparison','NumberTitle','off');
hold on

if has_lumped
    plot(x_stations, P0_L/1e3, 'o-', 'Color', col_L, 'LineWidth', lw, ...
         'MarkerFaceColor', col_L, 'MarkerSize', ms, 'DisplayName', 'Lumped')
end
plot(x_stations, P0_D/1e3, 's--', 'Color', col_D, 'LineWidth', lw, ...
     'MarkerFaceColor', col_D, 'MarkerSize', ms, 'DisplayName', 'Discretised')

set(gca, 'XTick', x_stations, 'XTickLabel', labels)
xlabel('TMS Station')
ylabel('Total pressure [kPa]')
title('Total pressure through TMS')
legend('Location', 'best')
grid on

%% ---- Figure 3: Static temperature through TMS ---- %%
figure('Name','Static Temperature — TMS Comparison','NumberTitle','off');
hold on

if has_lumped
    plot(x_stations, T_L, 'o-', 'Color', col_L, 'LineWidth', lw, ...
         'MarkerFaceColor', col_L, 'MarkerSize', ms, 'DisplayName', 'Lumped')
end
plot(x_stations, T_D, 's--', 'Color', col_D, 'LineWidth', lw, ...
     'MarkerFaceColor', col_D, 'MarkerSize', ms, 'DisplayName', 'Discretised')

% HX interior for discretised
plot(x_HX, T_air_seg, '-', 'Color', col_HX, 'LineWidth', lw, ...
     'DisplayName', 'Disc. HX interior')

set(gca, 'XTick', x_stations, 'XTickLabel', labels)
xlabel('TMS Station')
ylabel('Static temperature [K]')
title('Static temperature through TMS')
legend('Location', 'best')
grid on

%% ---- Figure 4: Velocity through TMS ---- %%
figure('Name','Velocity — TMS Comparison','NumberTitle','off');
hold on

if has_lumped
    plot(x_stations, v_L, 'o-', 'Color', col_L, 'LineWidth', lw, ...
         'MarkerFaceColor', col_L, 'MarkerSize', ms, 'DisplayName', 'Lumped')
end
plot(x_stations, v_D, 's--', 'Color', col_D, 'LineWidth', lw, ...
     'MarkerFaceColor', col_D, 'MarkerSize', ms, 'DisplayName', 'Discretised')

% HX interior for discretised
plot(x_HX, v_air_seg, '-', 'Color', col_HX, 'LineWidth', lw, ...
    'DisplayName', 'Disc. HX interior')
plot(x_HX, v_channel_seg, '--', 'Color', col_HX, 'LineWidth', lw, ...
    'DisplayName', 'Disc. HX interior channel')


set(gca, 'XTick', x_stations, 'XTickLabel', labels)
xlabel('TMS Station')
ylabel('Velocity [m/s]')
title('Velocity through TMS')
legend('Location', 'best')
grid on

%% ---- Figure 5: Mach number through TMS ---- %%
figure('Name','Mach Number — TMS Comparison','NumberTitle','off');
hold on

if has_lumped
    plot(x_stations, M_L, 'o-', 'Color', col_L, 'LineWidth', lw, ...
         'MarkerFaceColor', col_L, 'MarkerSize', ms, 'DisplayName', 'Lumped')
end
plot(x_stations, M_D, 's--', 'Color', col_D, 'LineWidth', lw, ...
     'MarkerFaceColor', col_D, 'MarkerSize', ms, 'DisplayName', 'Discretised')

set(gca, 'XTick', x_stations, 'XTickLabel', labels)
xlabel('TMS Station')
ylabel('Mach number [-]')
title('Mach number through TMS')
legend('Location', 'best')
grid on

%% ---- Figure 6: Combined 2x2 summary panel ---- %%
figure('Name','TMS Summary — Lumped vs Discretised','NumberTitle','off');

subplot(2,2,1)
hold on

if has_lumped
    plot(x_stations, P_L/1e3, 'o-', 'Color', col_L, 'LineWidth', lw, ...
         'MarkerFaceColor', col_L, 'MarkerSize', ms, 'DisplayName', 'Lumped')
end
plot(x_stations, P_D/1e3, 's--', 'Color', col_D, 'LineWidth', lw, ...
     'MarkerFaceColor', col_D, 'MarkerSize', ms, 'DisplayName', 'Discretised')
plot(x_HX, P_air_seg/1e3, '-', 'Color', col_HX, 'LineWidth', lw, ...
     'DisplayName', 'Disc. HX interior')
set(gca, 'XTick', x_stations, 'XTickLabel', labels)
ylabel('Static pressure [kPa]')
title('Static pressure')
legend('Location','best','FontSize',8)
grid on

subplot(2,2,2)
hold on

if has_lumped
    plot(x_stations, T_L, 'o-', 'Color', col_L, 'LineWidth', lw, ...
         'MarkerFaceColor', col_L, 'MarkerSize', ms, 'DisplayName', 'Lumped')
end
plot(x_stations, T_D, 's--', 'Color', col_D, 'LineWidth', lw, ...
     'MarkerFaceColor', col_D, 'MarkerSize', ms, 'DisplayName', 'Discretised')
plot(x_HX, T_air_seg, '-', 'Color', col_HX, 'LineWidth', lw, ...
     'DisplayName', 'Disc. HX interior')
set(gca, 'XTick', x_stations, 'XTickLabel', labels)
ylabel('Static temperature [K]')
title('Static temperature')
legend('Location','best','FontSize',8)
grid on

subplot(2,2,3)
hold on

if has_lumped
    plot(x_stations, v_L, 'o-', 'Color', col_L, 'LineWidth', lw, ...
         'MarkerFaceColor', col_L, 'MarkerSize', ms, 'DisplayName', 'Lumped')
end
plot(x_stations, v_D, 's--', 'Color', col_D, 'LineWidth', lw, ...
     'MarkerFaceColor', col_D, 'MarkerSize', ms, 'DisplayName', 'Discretised')
plot(x_HX, v_air_seg, '-', 'Color', col_HX, 'LineWidth', lw, ...
    'DisplayName', 'Disc. HX interior')
set(gca, 'XTick', x_stations, 'XTickLabel', labels)
ylabel('Velocity [m/s]')
title('Velocity')
legend('Location','best','FontSize',8)
grid on

subplot(2,2,4)
hold on

if has_lumped
    plot(x_stations, M_L, 'o-', 'Color', col_L, 'LineWidth', lw, ...
         'MarkerFaceColor', col_L, 'MarkerSize', ms, 'DisplayName', 'Lumped')
end
plot(x_stations, M_D, 's--', 'Color', col_D, 'LineWidth', lw, ...
     'MarkerFaceColor', col_D, 'MarkerSize', ms, 'DisplayName', 'Discretised')
set(gca, 'XTick', x_stations, 'XTickLabel', labels)
ylabel('Mach number [-]')
title('Mach number')
legend('Location','best','FontSize',8)
grid on

sgtitle('TMS Flow Variables — Lumped vs Discretised Model', 'FontSize', 12)

end