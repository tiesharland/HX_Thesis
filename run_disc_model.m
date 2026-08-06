function results = run_disc_model(N_segments, e, r, hx_theta, fan, fpr_init, ...
    M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
    d2_init, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, tol_T)

counter = 1;

%%%--- SHARED: FREESTREAM ---%%%

[P_inf, T_inf, Rho_inf, Gamma_inf, Cp_inf, mu_inf, P_inf_tot, M_inf, a_inf] = ...
    Freestream(h, p11, R, t11, V_inf);

P_inf_tot = P_inf + (0.5*Rho_inf*V_inf*V_inf);
a_inf     = sqrt(Gamma_inf*R*T_inf);
M_inf     = V_inf/a_inf;

%%%--- SHARED: PROPELLER ---%%%

[J_val, Cp_val, eff_val, D, J, Cp, V_ind_prop, dp_tot_prop, T_net, dia_prop, ...
    M1, P1, P1_0, T1, T1_0, Rho_1, v1, Re_1] = ...
    propeller(flight_phase, V_inf, Rho_inf, P_inf, a_inf, T_inf, mu_inf);

%%%--- SHARED: DIFFUSER ---%%%

A2             = pi*d2_init*d2_init/4;
A_frontal_target = A2*AR_diff;
d3_target      = sqrt(A_frontal_target*4/pi);
AR_diff_arr    = (d3_target/d2_init)^2;
AR_init        = AR_diff_arr;
theta_max_diff = 15;
flag           = "Friction";
d3             = d3_target;
m_dot_streamtube = Rho_1*v1*A2;
M_dot_in       = m_dot_streamtube;

[M3, P3, P3_0, A3, T3, T3_0, Rho_3, v3, M2, P2, P2_0, A2, T2, T2_0, Rho_2, v2, ...
    L_diffuser, M_dot_2, M_dot_3, d2, d3, Re_2, Re_3] = ...
    Diffuser_og_mine_for_7x7(M_dot_in, d3, M_dot_FOD, M_dot_comp, M1, P1, P1_0, ...
    T1, T1_0, Rho_1, v1, R, theta_max_diff, d2_init, A2, AR_init, a_inf, ...
    Gamma_inf, flag, mu_inf);

%%%--- SHARED: HX BOUNDARY CONDITIONS ---%%%

[T_mean_h, T_mean_c, T_c_i, T_c_o, T_h_i, T_h_o, C_h, C_c, C_star] = ...
    HX_deltaT(T3, M_dot_3, M_dot_coolant, T_out_fc, T_in_fc);

Pr_coolant = mu_EG(T_mean_h)*cp_EG_50_50(T_mean_h)/k_EG(T_mean_h);
Pr_air     = mu_air(T_mean_c)*cp_air(T_mean_c)/k_air(T_mean_c);

[dp_coolant_loop, d_h_air, M_dot_4, b_t_air, b_t_coolant, dp_hx, ...
          N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, NTU, R_tot, ...
          v_channel_air, v_channel_coolant, d_h_coolant, A_o_coolant, A_o_air, ...
          Re_air, Re_coolant, h_air, h_coolant, L_solution, UA_unit, v4, P4_0, ...
          M4, T4, T4_0, P4, F_drag_hx, M_hx, A4, drag_HX, T_cool_out_arr, dp_cool_seg, ...
          T_air_seg, P_air_seg, v_air_seg, Re_air_seg, Pr_air_seg, v_channel_seg, ...
          K_seg, f_air_seg, Nu_air_seg, h_air_seg, eta_fin_seg, ...
          UA_seg_arr, NTU_seg_arr, eps_seg_arr, Q_seg_arr, dp_seg_arr, ...
          Q_pred_solution, T_h_o_solution, P_0_air_seg, T_0_air_seg, M_air_seg, f_hx_seg, inlet_dp] = ...
          HX_design1_disc(e, r, hx_theta, counter, A3, v3, R, P3, ...
          d3, T_h_o, n_modules, T_h_i, T_c_i, T_mean_h, ...
          Q_tot, M_dot_coolant, M_dot_3, T3, N_segments, tol_T);

%%%--- PULLER FAN ---%%%

fpr_init = 1; % fan OFF

[M5, T5_0, T5, P5, P5_0, v5, A5, dp_fan, P_shaft, M_dot_5] = ...
    fan_backup(fpr_init, d3, P4, P4_0, T4, T4_0, v4, M4, A4, R, M_dot_4);

%%%--- NOZZLE ---%%%

% AR_noz = 0.33;
[A6, M6, L_nozzle, P6, P6_0, T6, T6_0, v6, M_dot_6, Re_5, d6, mdot_choked] = ...
    nozzle_old(AR_noz, flag, d3, P5, P5_0, T5, T5_0, v5, M5, M_dot_4);

%%%--- PRESSURE BALANCE ITERATION ---%%%

% With DeltaT
[P1, P2, P3, P4, P5, P6, P1_0, P2_0, P3_0, P4_0, P5_0, P6_0, ...
 M1, M2, M3, M4, M5, M6, v1, v2, v3, v4, v5, v6, ...
 T1, T1_0, T2, T2_0, T3, T3_0, T4, T4_0, T5, T5_0, T6, T6_0, ...
 M_dot_2, M_dot_3, M_dot_4, M_dot_5, M_dot_6, P_shaft, M_hx, dp_coolant, ...
 T_air_seg, P_air_seg, v_air_seg, Re_air_seg, Pr_air_seg, v_channel_seg, ...
 f_air_seg, Nu_air_seg, h_air_seg, Q_seg_arr, dp_seg_arr, T_cool_out_arr, dp_cool_seg, ...
 L_solution, Q_pred_solution, T_h_o_solution, dp_hx, drag_HX, ...
 P_0_air_seg, T_0_air_seg, M_air_seg, f_hx_seg, inlet_dp] = ...
    pressure_iteration_disc(e, r, m_dot_streamtube, fan, ...
    hx_theta, fpr_init, M_dot_in, flight_phase, theta_max_diff, d3, M_dot_FOD, ...
    M_dot_comp, counter, d2_init, AR_noz, P1, P_inf, ...
    T_inf, V_inf, Rho_inf, AR_init, a_inf, Gamma_inf, flag, mu_inf, ...
    A2, v2, M2, R, P2, P2_0, T2_0, n_modules, Q_tot, M_dot_coolant, ...
    M_dot_2, M_dot_3, T2, P3, P3_0, T3, T3_0, v3, M3, P4, P4_0, T4, T4_0, ...
    v4, M4, P5, P5_0, M5, T5, T5_0, v5, P6, P6_0, M6, T6, T6_0, v6, ...
    M1, T1, T1_0, v1, P1_0, T_out_fc, T_in_fc, N_segments, tol_T);


if fpr_init == 1
    fprintf("The puller fan is not present\n")
else 
    fprintf("The puller fan is present\n")
end


m_dot_seg = zeros(1, N_segments+1);
for k = 1:N_segments+1
    rho_k        = rho_air(P_air_seg(k), T_air_seg(k));
    m_dot_seg(k) = rho_k * v_channel_seg(k) * A_o_air * N_fin_air * N_air_pass;
end
m_dot_seg(1) = m_dot_seg(2);


m_dot_spill = m_dot_streamtube-M_dot_2;
fprintf("Spillage mass flow rate = %f\n", m_dot_spill);
drag_tot = (m_dot_streamtube*v1) - (M_dot_6*v6) - (M_dot_2*M_dot_FOD*v3) + ((P2-P_inf)*A2) - ((P6-P_inf)*A6);


fprintf("Total drag = %f\n",drag_tot);


results.L_solution   = L_solution;
results.dp_hx        = dp_hx;
results.dp_coolant   = dp_coolant;
results.T4           = T4;
results.M_hx         = M_hx;
results.drag_HX      = drag_HX;
results.drag_tot     = drag_tot;
results.T_cool_out   = T_h_o_solution;
results.P_inf        = P_inf;
results.P1           = P1;
results.P2           = P2;
results.P3           = P3;
results.P4           = P4;
results.P5           = P5;
results.P6           = P6;
results.T_inf        = T_inf;
results.T1           = T1;
results.T2           = T2;
results.T3           = T3;
results.T5           = T5;
results.T6           = T6;
results.P_inf_tot    = P_inf_tot;
results.P1_0         = P1_0;
results.P2_0         = P2_0;
results.P3_0         = P3_0;
results.P4_0         = P4_0;
results.P5_0         = P5_0;
results.P6_0         = P6_0;
results.V_inf        = V_inf;
results.v1           = v1;
results.v2           = v2;
results.v3           = v3;
results.v4           = v4;
results.v5           = v5;
results.v6           = v6;
results.M_inf        = M_inf;
results.M1           = M1;
results.M2           = M2;
results.M3           = M3;
results.M4           = M4;
results.M5           = M5;
results.M6           = M6;
results.M_dot_2      = M_dot_2;
results.M_dot_3      = M_dot_3;
results.M_dot_4      = M_dot_4;
results.M_dot_5      = M_dot_5;
results.M_dot_6      = M_dot_6;
results.m_spill      = m_dot_streamtube - M_dot_2;
results.m_dot_seg    = m_dot_seg;
results.AR_init      = AR_init;
results.AR_noz       = AR_noz;
results.fpr          = fpr_init;
results.T_c_i        = T_c_i;
results.T_h_i        = T_h_i;
results.T_h_o        = T_h_o;
results.N_segments   = N_segments;
results.P_air_seg    = P_air_seg;
results.T_air_seg    = T_air_seg;
results.Q_seg_arr    = Q_seg_arr;
results.dp_seg_arr   = dp_seg_arr;
results.v_air_seg    = v_air_seg;
results.v_channel_seg = v_channel_seg;
results.Re_air_seg   = Re_air_seg;
results.Pr_air_seg   = Pr_air_seg;
results.f_air_seg    = f_air_seg;
results.Nu_air_seg   = Nu_air_seg;
results.h_air_seg    = h_air_seg;
results.Q_pred_solution = Q_pred_solution;
results.T_cool_out_arr = T_cool_out_arr;
results.dp_cool_seg  = dp_cool_seg;
results.P_0_air_seg  = P_0_air_seg;
results.T_0_air_seg  = T_0_air_seg;
results.M_air_seg    = M_air_seg;
results.f_hx_seg     = f_hx_seg;
results.inlet_dp     = inlet_dp;

end