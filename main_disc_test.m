clear all
clc

counter = 1;

d2_init = 0.45;
AR_diff = 4;
A2      = pi*d2_init*d2_init/4;
A_frontal_target = A2*AR_diff;
d3_target        = sqrt(A_frontal_target*4/pi);
AR_diff_arr      = (d3_target/d2_init)^2;
n_modules        = 2;

%%%--- INPUTS - FREESTREAM ---%%%

flight_phase = 3; % 1=take-off, 2=top of climb, 3=cruise

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

p11    = 22632;
t11    = 216.65;
R      = 287;
hx_theta = 60;

fan  = "OFF";
disc = "ON";
N_segments = 10;

%%%--- MASS FLOW LEAKS ---%%%

M_dot_FOD  = 0;
M_dot_comp = 0;

%%%--- DIFFUSER INPUTS ---%%%

AR_init = AR_diff_arr;

%%%--- FREESTREAM ---%%%

[P_inf, T_inf, Rho_inf, Gamma_inf, Cp_inf, mu_inf, P_inf_tot, M_inf, a_inf] = ...
    Freestream(h, p11, R, t11, V_inf);

P_inf_tot = P_inf + (0.5*Rho_inf*V_inf*V_inf);
a_inf     = sqrt(Gamma_inf*R*T_inf);
M_inf     = V_inf/a_inf;

%%%--- PROPELLER ---%%%

[J_val, Cp_val, eff_val, D, J, Cp, V_ind_prop, dp_tot_prop, T_net, dia_prop, ...
 M1, P1, P1_0, T1, T1_0, Rho_1, v1, Re_1] = ...
    propeller(flight_phase, V_inf, Rho_inf, P_inf, a_inf, T_inf, mu_inf);

%%%--- DIFFUSER ---%%%

theta_max_diff   = 15;
flag             = "Friction";
d3               = d3_target;
m_dot_streamtube = Rho_1*v1*A2;
M_dot_in         = m_dot_streamtube;

[M3, P3, P3_0, A3, T3, T3_0, Rho_3, v3, M2, P2, P2_0, A2, T2, T2_0, Rho_2, v2, ...
 L_diffuser, M_dot_2, M_dot_3, d2, d3, Re_2, Re_3] = ...
    Diffuser_og_mine_for_7x7(M_dot_in, d3, M_dot_FOD, M_dot_comp, M1, P1, P1_0, ...
    T1, T1_0, Rho_1, v1, R, theta_max_diff, d2_init, A2, AR_init, a_inf, ...
    Gamma_inf, flag, mu_inf);

%%%--- HX INPUTS ---%%%

M_dot_coolant_max = 44.4;
M_dot_coolant_min = 29;

M_dot_coolant = 44.4;

%%%--- MEAN TEMPERATURE CALCULATIONS ---%%%

[T_mean_h, T_mean_c, T_c_i, T_c_o, T_h_i, T_h_o, C_h, C_c, C_star] = ...
    HX_deltaT(T3, M_dot_3, M_dot_coolant, T_out_fc, T_in_fc);

%%%--- PRANDTL NUMBERS ---%%%

Pr_coolant = mu_EG(T_mean_h)*cp_EG_50_50(T_mean_h)/k_EG(T_mean_h);

%%%--- HX SIZING ---%%%

e = 4.4;
r = 3.8;
tol_T = 7.5; % 8.0557;

[dp_coolant_loop, d_h_air, M_dot_4, b_t_air, b_t_coolant, dp_hx, ...
    N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, NTU, R_tot, ...
    v_channel_air, v_channel_coolant, d_h_coolant, A_o_coolant, A_o_air, ...
    Re_air, Re_coolant, h_air, h_coolant, L_solution, UA_unit, v4, P4_0, ...
    M4, T4, T4_0, P4, F_drag_hx, M_hx, A4, drag_HX, T_cool_out_arr, dp_cool_seg, ...
    T_air_seg, P_air_seg, v_air_seg, Re_air_seg, Pr_air_seg, v_channel_seg, ...
    K_seg, f_air_seg, Nu_air_seg, h_air_seg, eta_fin_seg, ...
    UA_seg_arr, NTU_seg_arr, eps_seg_arr, Q_seg_arr, dp_seg_arr, ...
    Q_pred_solution, T_h_o_solution] = ...
HX_design1_disc(e, r, hx_theta, counter, A3, v3, R, P3, ...
    d3, T_h_o, n_modules, T_h_i, T_c_i, T_mean_h, ...
    Q_tot, M_dot_coolant, M_dot_3, T3, N_segments, tol_T);

%%%--- PULLER FAN ---%%%

fpr_init = 1; % fan OFF

[M5, T5_0, T5, P5, P5_0, v5, A5, dp_fan, P_shaft, M_dot_5] = ...
    fan_backup(fpr_init, d3, P4, P4_0, T4, T4_0, v4, M4, A4, R, M_dot_4);

%%%--- NOZZLE ---%%%

AR_noz = 0.33;
[A6, M6, L_nozzle, P6, P6_0, T6, T6_0, v6, M_dot_6, Re_5, d6, mdot_choked] = ...
    nozzle_old(AR_noz, flag, d3, P5, P5_0, T5, T5_0, v5, M5, M_dot_4);

%%%--- PRESSURE BALANCE ITERATION ---%%%

% With DeltaT
[P1, P2, P3, P4, P5, P6, P1_0, P2_0, P3_0, P4_0, P5_0, P6_0, ...
    M1, M2, M3, M4, M5, M6, v1, v2, v3, v4, v5, v6, ...
    T1, T1_0, T2, T2_0, T3, T3_0, T4, T4_0, T5, T5_0, T6, T6_0, ...
    M_dot_2, M_dot_3, M_dot_4, M_dot_6, P_shaft, M_hx, dp_coolant, ...
    T_air_seg, P_air_seg, v_air_seg, Re_air_seg, v_channel_seg, ...
    f_air_seg, Nu_air_seg, h_air_seg, Q_seg_arr, dp_seg_arr, T_cool_out_arr, dp_cool_seg, ...
    L_solution, Q_pred_solution, T_h_o_solution, dp_hx, drag_HX] = ...
pressure_iteration_disc(e, r, m_dot_streamtube, fan, ...
    hx_theta, fpr_init, M_dot_in, flight_phase, theta_max_diff, d3, M_dot_FOD, ...
    M_dot_comp, counter, d2_init, AR_noz, P1, P_inf, ...
    T_inf, V_inf, Rho_inf, AR_init, a_inf, Gamma_inf, flag, mu_inf, ...
    A2, v2, M2, R, P2, P2_0, T2_0, n_modules, Q_tot, M_dot_coolant, ...
    M_dot_2, M_dot_3, T2, P3, P3_0, T3, T3_0, v3, M3, P4, P4_0, T4, T4_0, ...
    v4, M4, P5, P5_0, M5, T5, T5_0, v5, P6, P6_0, M6, T6, T6_0, v6, ...
    M1, T1, T1_0, v1, P1_0, T_out_fc, T_in_fc, N_segments, tol_T);

% % No DeltaT
% [P1, P2, P3, P4, P5, P6, P1_0, P2_0, P3_0, P4_0, P5_0, P6_0, ...
% M1, M2, M3, M4, M5, M6, v1, v2, v3, v4, v5, v6, ...
% T1, T1_0, T2, T2_0, T3, T3_0, T4, T4_0, T5, T5_0, T6, T6_0, ...
% M_dot_2, M_dot_3, M_dot_4, M_dot_6, P_shaft, M_hx, dp_coolant_loop, ...
% T_air_seg, P_air_seg, v_air_seg, Re_air_seg, ...
% f_air_seg, Nu_air_seg, h_air_seg, Q_seg_arr, dp_seg_arr] = ...
% pressure_iteration_disc_no_deltaT(e, r, m_dot_streamtube, fan, hx_theta, fpr_init, M_dot_in, flight_phase, theta_max_diff, ...
% d3, M_dot_FOD, M_dot_comp, counter, d2_init, AR_noz, P1, P_inf, T_inf, V_inf, Rho_inf, AR_init, ...
% a_inf, Gamma_inf, flag, mu_inf, A2, v2, M2, R, P2, P2_0, ...
% T2_0, n_modules, Q_tot, ...
% M_dot_coolant, M_dot_2, M_dot_3, T2, P3, ...
% P3_0, T3, T3_0, v3, M3, P4, P4_0, T4, T4_0, v4, M4, P5, P5_0, ...
% M5, T5, T5_0, v5, P6, P6_0, M6, T6, T6_0, v6, M1, T1, T1_0, ...
% v1, P1_0, T_out_fc, T_in_fc, N_segments, T_h_o, T_h_i, T_c_i, T_mean_h);

m_dot_seg = zeros(1, N_segments+1);
for k = 1:N_segments+1
    rho_k        = rho_air(P_air_seg(k), T_air_seg(k));
    m_dot_seg(k) = rho_k * v_channel_seg(k) * A_o_air * N_fin_air * N_air_pass;
end
m_dot_seg(1) = m_dot_seg(2);

%%%--- SUMMARY PRINTS ---%%%

fprintf("\n=== HX RESULTS ===\n");
fprintf("HX length          = %.4f m\n",   L_solution);
fprintf("HX pressure drop   = %.2f Pa\n",  dp_hx);
fprintf("Air outlet temp    = %.2f K\n",   T4);
fprintf("HX mass            = %.2f kg\n",  M_hx);
fprintf("HX drag         = %.2f N\n",   drag_HX);
fprintf("P6 - P_inf         = %.2f Pa\n",  P6 - P_inf);
fprintf("Coolant outlet temp = %.2f K\n",   T_h_o_solution);
fprintf("Coolant pressure drop (avg) = %.2f Pa\n",   dp_coolant);
fprintf("Cooling power  = %.2f MW\n",   Q_pred_solution/1e6);
fprintf("Nozzle exit P6     = %.2f Pa\n",  P6);
fprintf("Ambient P_inf      = %.2f Pa\n",  P_inf);

%%%--- PLOTS ---%%%

% plot1(AR_init, AR_noz, fpr_init, P_inf, P1, P2, P3, P4, P5,P6, T_inf,T1, T2, T3, T4, T5, T6, P_inf_tot, P1_0, P2_0, P3_0, P4_0, P5_0, P6_0, V_inf, v1, v2, v3, v4, v5, v6, M_inf, M1, M2, M3, M4, M5, M6);

plot_stations(P_inf, P_inf_tot, T_inf, M_inf, V_inf, ...
    P1, P2, P3, P4, P5, P6, ...
    T1, T2, T3, T4, T5, T6, ...
    P1_0, P2_0, P3_0, P4_0, P5_0, P6_0, ...
    v1, v2, v3, v4, v5, v6, ...
    M1, M2, M3, M4, M5, M6, ...
    N_segments, T_air_seg, P_air_seg, v_air_seg, v_channel_seg);


if fpr_init == 1
    fprintf("The puller fan is not present\n")
else 
    fprintf("The puller fan is present\n")
end


m_dot_spill = m_dot_streamtube-M_dot_2;
fprintf("Spillage mass flow rate = %f\n", m_dot_spill);
drag_tot = (m_dot_streamtube*v1) - (M_dot_6*v6) - (M_dot_2*M_dot_FOD*v3) + ((P2-P_inf)*A2) - ((P6-P_inf)*A6);


fprintf("Total drag = %f\n",drag_tot);

plot_HX(N_segments, L_solution, T_h_i, T_c_i, ...
    T_air_seg, P_air_seg, v_air_seg, Re_air_seg, ...
    f_air_seg, Nu_air_seg, h_air_seg, Q_seg_arr, dp_seg_arr, ...
    T_cool_out_arr, dp_cool_seg, ...
    P3, P4, M_dot_3, M_dot_4, dp_hx, m_dot_seg);