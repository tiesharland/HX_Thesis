clear all
clc

%%%--- SHARED INPUTS ---%%%

flight_phase = 3;

if flight_phase == 1
    h        = 5;
    T_in_fc  = 85+273;
    T_out_fc = 105+273;
    Q_tot    = 2.38e6;
    V_inf    = 2;
elseif flight_phase == 3
    h        = 4876;
    T_in_fc  = 70+273;
    T_out_fc = 85+273;
    V_inf    = 128;
    Q_tot    = 2.25e6;
elseif flight_phase == 2
    h        = 4876;
    T_in_fc  = 70+273;
    T_out_fc = 85+273;
    Q_tot    = 2.25e6;
    V_inf    = 128;
end

p11           = 22632;
t11           = 216.65;
R             = 287;
hx_theta      = 60;
fan           = "OFF";
e             = 4.4;
r             = 3.8;
% n_modules     = 2;
d2_init       = 0.45;
AR_diff       = 4;
AR_noz        = 0.33;
fpr_init      = 1;
M_dot_coolant = 44.4; % Max: 44.4, Min: 29
N_segments    = 10;

M_dot_FOD  = 0;
M_dot_comp = 0;

% d2_init = 0.45; %m, diffuser inlet initial guess
A2 = pi*d2_init*d2_init/4;
A_frontal_target = A2*AR_diff; %this number multiplied to A2 is your diffuser AR
d3_target = sqrt(A_frontal_target.*4/pi);
AR_diff_arr = (d3_target./d2_init).^2;
n_modules = 2; %ducts per nacelle

%%%--- RUN LUMPED MODEL ---%%%

fprintf("\n========================================\n");
fprintf("        RUNNING LUMPED MODEL            \n");
fprintf("========================================\n");

results_L = run_lumped_model(e, r, hx_theta, fan, fpr_init, ...
    M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
    d2_init, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, AR_diff_arr);

%%%--- RUN DISCRETISED MODEL ---%%%

fprintf("\n========================================\n");
fprintf("      RUNNING DISCRETISED MODEL         \n");
fprintf("========================================\n");

% tol_T = results_L.T_cool_out - results_L.T_h_o;
tol_T = 7.5;

results_D = run_disc_model(N_segments, e, r, hx_theta, fan, fpr_init, ...
    M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11,...
    d2_init, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, tol_T);

%%%--- COMPARISON SUMMARY ---%%%

fprintf("\n========================================\n");
fprintf("           MODEL COMPARISON             \n");
fprintf("========================================\n");
fprintf("%-30s %-15s %-15s\n", "Quantity", "Lumped", "Discretised");
fprintf("%-30s %-15.4f %-15.4f\n", "HX length [m]",       results_L.L_solution,  results_D.L_solution);
fprintf("%-30s %-15.2f %-15.2f\n", "HX dp [Pa]",          results_L.dp_hx,       results_D.dp_hx);
fprintf("%-30s %-15.2f %-15.2f\n", "Air outlet T4 [K]",   results_L.T4,          results_D.T4);
fprintf("%-30s %-15.2f %-15.2f\n", "HX mass [kg]",        results_L.M_hx,        results_D.M_hx);
fprintf("%-30s %-15.2f %-15.2f\n", "HX drag [N]",         results_L.drag_HX,     results_D.drag_HX);
fprintf("%-30s %-15.2f %-15.2f\n", "P6 - P_inf [Pa]",     results_L.P6-results_L.P_inf,    results_D.P6-results_D.P_inf);
fprintf("%-30s %-15.2f %-15.2f\n", "Coolant outlet temp [K]",  results_L.T_cool_out,  results_D.T_cool_out);
fprintf("%-30s %-15.2f %-15.2f\n", "Coolant pressure drop [K]",  results_L.dp_coolant,  results_D.dp_coolant);
fprintf("%-30s %-15.4f %-15.4f\n", "Cooling power [MW]", results_L.Q_pred_solution/1e6,     results_D.Q_pred_solution/1e6);
fprintf("%-30s %-15.4f %-15.4f\n", "Mass flow in [kg/s]", results_L.M_dot_2,     results_D.M_dot_2);
fprintf("%-30s %-15.4f %-15.4f\n", "Spillage [kg/s]",     results_L.m_spill,     results_D.m_spill);

%%%--- PLOTS ---%%%

plot_stations(results_D.P_inf, results_D.P_inf_tot, results_D.T_inf, results_D.M_inf, results_D.V_inf, ...
    results_D.P1, results_D.P2, results_D.P3, results_D.P4, results_D.P5, results_D.P6, ...
    results_D.T1, results_D.T2, results_D.T3, results_D.T4, results_D.T5, results_D.T6, ...
    results_D.P1_0, results_D.P2_0, results_D.P3_0, results_D.P4_0, results_D.P5_0, results_D.P6_0, ...
    results_D.v1, results_D.v2, results_D.v3, results_D.v4, results_D.v5, results_D.v6, ...
    results_D.M1, results_D.M2, results_D.M3, results_D.M4, results_D.M5, results_D.M6, ...
    N_segments, results_D.T_air_seg, results_D.P_air_seg, results_D.v_air_seg, results_D.v_channel_seg, ...
    results_L.P1, results_L.P2, results_L.P3, results_L.P4, results_L.P5, results_L.P6, ...
    results_L.T1, results_L.T2, results_L.T3, results_L.T4, results_L.T5, results_L.T6, ...
    results_L.P1_0, results_L.P2_0, results_L.P3_0, results_L.P4_0, results_L.P5_0, results_L.P6_0, ...
    results_L.v1, results_L.v2, results_L.v3, results_L.v4, results_L.v5, results_L.v6, ...
    results_L.M1, results_L.M2, results_L.M3, results_L.M4, results_L.M5, results_L.M6);

plot_HX(N_segments, results_D.L_solution, results_D.T_h_i, results_D.T_c_i, ...
    results_D.T_air_seg, results_D.P_air_seg, results_D.v_air_seg, results_D.Re_air_seg, ...
    results_D.f_air_seg, results_D.Nu_air_seg, results_D.h_air_seg, results_D.Q_seg_arr, results_D.dp_seg_arr, ...
    results_D.T_cool_out_arr, results_D.dp_cool_seg, ...
    results_D.P3, results_D.P4, results_D.M_dot_3, results_D.M_dot_4, results_D.dp_hx, results_D.m_dot_seg, ...
    results_L.L_solution, results_L.T4, results_L.P3, results_L.P4, results_L.T3, results_L.dp_hx);
