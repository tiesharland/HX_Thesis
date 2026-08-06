function plot_stations(results_D, varargin)
% Usage:
%   plot_stations(results_D)                          — discretised only
%   plot_stations(results_D, results_L)                — + lumped comparison
%   plot_stations(results_D, results_L, results_DNS)   — + DNS comparison
%   plot_stations(results_D, [], results_DNS)          — DNS only, no lumped

has_lumped = ~isempty(varargin) && ~isempty(varargin{1});
has_dns    = numel(varargin) >= 2 && ~isempty(varargin{2});

if has_lumped
    results_L = varargin{1};
end
if has_dns
    results_DNS = varargin{2};
end

%% ---- Shared freestream quantities (from results_D only) ---- %%

P_inf     = results_D.P_inf;
P_inf_tot = results_D.P_inf_tot;
T_inf     = results_D.T_inf;
M_inf     = results_D.M_inf;
V_inf     = results_D.V_inf;

%% ---- TMS state arrays — discretised ---- %%

P_D  = [P_inf,     results_D.P1,   results_D.P2,   results_D.P3,   results_D.P4,   results_D.P5,   results_D.P6];
P0_D = [P_inf_tot, results_D.P1_0, results_D.P2_0, results_D.P3_0, results_D.P4_0, results_D.P5_0, results_D.P6_0];
T_D  = [T_inf,     results_D.T1,   results_D.T2,   results_D.T3,   results_D.T4,   results_D.T5,   results_D.T6];
v_D  = [V_inf,     results_D.v1,   results_D.v2,   results_D.v3,   results_D.v4,   results_D.v5,   results_D.v6];
M_D  = [M_inf,     results_D.M1,   results_D.M2,   results_D.M3,   results_D.M4,   results_D.M5,   results_D.M6];

N_segments    = results_D.N_segments;
T_air_seg     = results_D.T_air_seg;
P_air_seg     = results_D.P_air_seg;
P_0_air_seg   = results_D.P_0_air_seg;
v_air_seg     = results_D.v_air_seg;
v_channel_seg = results_D.v_channel_seg;

%% ---- TMS state arrays — lumped ---- %%

if has_lumped
    P_L  = [P_inf,     results_L.P1,   results_L.P2,   results_L.P3,   results_L.P4,   results_L.P5,   results_L.P6];
    P0_L = [P_inf_tot, results_L.P1_0, results_L.P2_0, results_L.P3_0, results_L.P4_0, results_L.P5_0, results_L.P6_0];
    T_L  = [T_inf,     results_L.T1,   results_L.T2,   results_L.T3,   results_L.T4,   results_L.T5,   results_L.T6];
    v_L  = [V_inf,     results_L.v1,   results_L.v2,   results_L.v3,   results_L.v4,   results_L.v5,   results_L.v6];
    M_L  = [M_inf,     results_L.M1,   results_L.M2,   results_L.M3,   results_L.M4,   results_L.M5,   results_L.M6];
end

%% ---- TMS state arrays — DNS-based ---- %%

if has_dns
    P_DNS  = [P_inf,     results_DNS.P1,   results_DNS.P2,   results_DNS.P3,   results_DNS.P4,   results_DNS.P5,   results_DNS.P6];
    P0_DNS = [P_inf_tot, results_DNS.P1_0, results_DNS.P2_0, results_DNS.P3_0, results_DNS.P4_0, results_DNS.P5_0, results_DNS.P6_0];
    T_DNS  = [T_inf,     results_DNS.T1,   results_DNS.T2,   results_DNS.T3,   results_DNS.T4,   results_DNS.T5,   results_DNS.T6];
    v_DNS  = [V_inf,     results_DNS.v1,   results_DNS.v2,   results_DNS.v3,   results_DNS.v4,   results_DNS.v5,   results_DNS.v6];
    M_DNS  = [M_inf,     results_DNS.M1,   results_DNS.M2,   results_DNS.M3,   results_DNS.M4,   results_DNS.M5,   results_DNS.M6];

    N_segments_DNS    = results_DNS.N_segments;
    T_air_seg_DNS     = results_DNS.T_air_seg;
    P_air_seg_DNS     = results_DNS.P_air_seg;
    P_0_air_seg_DNS   = results_DNS.P_0_air_seg;
    v_air_seg_DNS     = results_DNS.v_air_seg;
    v_channel_seg_DNS = results_DNS.v_channel_seg;
end

%% ---- Colour scheme ---- %%
col_L   = [0.122 0.471 0.706];   % blue   — lumped
col_D   = [0.839 0.153 0.157];   % red    — discretised
col_HX  = [0.173 0.627 0.173];   % green  — HX interior (discretised)
col_DNS = [0.580 0.404 0.741];   % purple — DNS-based
col_HXDNS = [1.000 0.498 0.055]; % orange — HX interior (DNS-based)
lw = 1.5;
ms = 7;

%% ---- Station x-positions ---- %%
x_stations = 0:6;
labels     = {'0 (∞)','1','2','3','4','5','6'};

x_HX     = linspace(3, 4, N_segments+1);
if has_dns
    x_HX_DNS = linspace(3, 4, N_segments_DNS+1);
end

%% ---- Figure 1: Static pressure through TMS ---- %%
figure('Name','Static Pressure — TMS Comparison','NumberTitle','off');
hold on

if has_lumped
    plot(x_stations, P_L/1e3, 'o-', 'Color', col_L, 'LineWidth', lw, ...
         'MarkerFaceColor', col_L, 'MarkerSize', ms, 'DisplayName', 'Lumped')
end
plot(x_stations, P_D/1e3, 's--', 'Color', col_D, 'LineWidth', lw, ...
     'MarkerFaceColor', col_D, 'MarkerSize', ms, 'DisplayName', 'Discretised')
plot(x_HX, P_air_seg/1e3, '-', 'Color', col_HX, 'LineWidth', lw, ...
     'DisplayName', 'Disc. HX interior')
if has_dns
    plot(x_stations, P_DNS/1e3, '^:', 'Color', col_DNS, 'LineWidth', lw, ...
         'MarkerFaceColor', col_DNS, 'MarkerSize', ms, 'DisplayName', 'DNS-based')
    plot(x_HX_DNS, P_air_seg_DNS/1e3, '-', 'Color', col_HXDNS, 'LineWidth', lw, ...
         'DisplayName', 'DNS HX interior')
end

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
plot(x_HX, P_0_air_seg/1e3, '-', 'Color', col_HX, 'LineWidth', lw, ...
     'DisplayName', 'Disc. HX interior')
if has_dns
    plot(x_stations, P0_DNS/1e3, '^:', 'Color', col_DNS, 'LineWidth', lw, ...
         'MarkerFaceColor', col_DNS, 'MarkerSize', ms, 'DisplayName', 'DNS-based')
    plot(x_HX_DNS, P_0_air_seg_DNS/1e3, '-', 'Color', col_HXDNS, 'LineWidth', lw, ...
         'DisplayName', 'DNS HX interior')
end

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
plot(x_HX, T_air_seg, '-', 'Color', col_HX, 'LineWidth', lw, ...
     'DisplayName', 'Disc. HX interior')
if has_dns
    plot(x_stations, T_DNS, '^:', 'Color', col_DNS, 'LineWidth', lw, ...
         'MarkerFaceColor', col_DNS, 'MarkerSize', ms, 'DisplayName', 'DNS-based')
    plot(x_HX_DNS, T_air_seg_DNS, '-', 'Color', col_HXDNS, 'LineWidth', lw, ...
         'DisplayName', 'DNS HX interior')
end

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
plot(x_HX, v_air_seg, '-', 'Color', col_HX, 'LineWidth', lw, ...
     'DisplayName', 'Disc. HX interior')
plot(x_HX, v_channel_seg, '--', 'Color', col_HX, 'LineWidth', lw, ...
     'DisplayName', 'Disc. HX interior channel')
if has_dns
    plot(x_stations, v_DNS, '^:', 'Color', col_DNS, 'LineWidth', lw, ...
         'MarkerFaceColor', col_DNS, 'MarkerSize', ms, 'DisplayName', 'DNS-based')
    plot(x_HX_DNS, v_air_seg_DNS, '-', 'Color', col_HXDNS, 'LineWidth', lw, ...
         'DisplayName', 'DNS HX interior')
    plot(x_HX_DNS, v_channel_seg_DNS, '--', 'Color', col_HXDNS, 'LineWidth', lw, ...
         'DisplayName', 'DNS HX interior channel')
end

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
if has_dns
    plot(x_stations, M_DNS, '^:', 'Color', col_DNS, 'LineWidth', lw, ...
         'MarkerFaceColor', col_DNS, 'MarkerSize', ms, 'DisplayName', 'DNS-based')
end

set(gca, 'XTick', x_stations, 'XTickLabel', labels)
xlabel('TMS Station')
ylabel('Mach number [-]')
title('Mach number through TMS')
legend('Location', 'best')
grid on

%% ---- Figure 6: Combined 2x2 summary panel ---- %%
figure('Name','TMS Summary — Lumped vs Discretised vs DNS','NumberTitle','off');

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
if has_dns
    plot(x_stations, P_DNS/1e3, '^:', 'Color', col_DNS, 'LineWidth', lw, ...
         'MarkerFaceColor', col_DNS, 'MarkerSize', ms, 'DisplayName', 'DNS-based')
    plot(x_HX_DNS, P_air_seg_DNS/1e3, '-', 'Color', col_HXDNS, 'LineWidth', lw, ...
         'DisplayName', 'DNS HX interior')
end
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
if has_dns
    plot(x_stations, T_DNS, '^:', 'Color', col_DNS, 'LineWidth', lw, ...
         'MarkerFaceColor', col_DNS, 'MarkerSize', ms, 'DisplayName', 'DNS-based')
    plot(x_HX_DNS, T_air_seg_DNS, '-', 'Color', col_HXDNS, 'LineWidth', lw, ...
         'DisplayName', 'DNS HX interior')
end
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
if has_dns
    plot(x_stations, v_DNS, '^:', 'Color', col_DNS, 'LineWidth', lw, ...
         'MarkerFaceColor', col_DNS, 'MarkerSize', ms, 'DisplayName', 'DNS-based')
    plot(x_HX_DNS, v_air_seg_DNS, '-', 'Color', col_HXDNS, 'LineWidth', lw, ...
         'DisplayName', 'DNS HX interior')
end
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
if has_dns
    plot(x_stations, M_DNS, '^:', 'Color', col_DNS, 'LineWidth', lw, ...
         'MarkerFaceColor', col_DNS, 'MarkerSize', ms, 'DisplayName', 'DNS-based')
end
set(gca, 'XTick', x_stations, 'XTickLabel', labels)
ylabel('Mach number [-]')
title('Mach number')
legend('Location','best','FontSize',8)
grid on

sgtitle('TMS Flow Variables — Lumped vs Discretised vs DNS Model', 'FontSize', 12)

%% ---- Figure 7: Mass flow through TMS ---- %%
figure('Name','Mass Flow — TMS Comparison','NumberTitle','off');
hold on

% Reconstruct station mass flows from available fields
m_dot_streamtube_D = results_D.M_dot_2 + results_D.m_spill;
M_dot_D = [m_dot_streamtube_D, m_dot_streamtube_D, ...
    results_D.M_dot_2, results_D.M_dot_3, results_D.M_dot_4, ...
    results_D.M_dot_5, results_D.M_dot_6];

if has_lumped
    m_dot_streamtube_L = results_L.M_dot_2 + results_L.m_spill;
    M_dot_L = [m_dot_streamtube_L, m_dot_streamtube_L, ...
        results_L.M_dot_2, results_L.M_dot_3, results_L.M_dot_4, ...
        results_L.M_dot_5, results_L.M_dot_6];
end

if has_dns
    m_dot_streamtube_DNS = results_DNS.M_dot_2 + results_DNS.m_spill;
    M_dot_DNS = [m_dot_streamtube_DNS, m_dot_streamtube_DNS, ...
        results_DNS.M_dot_2, results_DNS.M_dot_3, results_DNS.M_dot_4, ...
        results_DNS.M_dot_5, results_DNS.M_dot_6];
end

if has_lumped
    plot(x_stations, M_dot_L, 'o-', 'Color', col_L, 'LineWidth', lw, ...
        'MarkerFaceColor', col_L, 'MarkerSize', ms, 'DisplayName', 'Lumped')
end

plot(x_stations, M_dot_D, 's--', 'Color', col_D, 'LineWidth', lw, ...
    'MarkerFaceColor', col_D, 'MarkerSize', ms, 'DisplayName', 'Discretised')

% HX interior mass flow through core nodes
plot(x_HX, results_D.m_dot_seg, '-', 'Color', col_HX, 'LineWidth', lw, ...
    'DisplayName', 'Disc. HX interior')

if has_dns
    plot(x_stations, M_dot_DNS, '^:', 'Color', col_DNS, 'LineWidth', lw, ...
        'MarkerFaceColor', col_DNS, 'MarkerSize', ms, 'DisplayName', 'DNS-based')
    plot(x_HX_DNS, results_DNS.m_dot_seg, '-', 'Color', col_HXDNS, 'LineWidth', lw, ...
        'DisplayName', 'DNS HX interior')
end

set(gca, 'XTick', x_stations, 'XTickLabel', labels)
xlabel('TMS Station')
ylabel('Mass flow rate [kg/s]')
title('Mass flow rate through TMS')
legend('Location', 'best')
grid on

end