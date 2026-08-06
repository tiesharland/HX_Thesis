function [P1,P2,P3,P4,P5,P6,P1_0,P2_0,P3_0,P4_0,P5_0,P6_0, M1, M2, M3, M4, M5, M6,  v1,v2,v3,v4,v5,v6, T1, T1_0, T2, T2_0, T3,T3_0,T4, T4_0, T5, T5_0, T6, T6_0, M_dot_2, M_dot_3, M_dot_4, M_dot_5, M_dot_6,P_shaft, M_hx, dp_coolant, L_solution, Q_pred_sol, T_cool_out, dp_hx, drag_HX, h_air, v_channel_air, Re_air, Pr_air, f_air, Nu_air, f_hx, inlet_dp] = pressure_eq_mDot_adjust_9thJan_noExhaust(e,r,m_dot_streamtube, fan, M_dot_des, FPR_des, hx_theta,fpr_init, M_dot_in, flight_phase, theta_max_diff,d3,M_dot_FOD, M_dot_comp,counter,A2_init, d2_init,AR_noz, dp_tot_prop, P1,P_inf, P_inf_tot, T_inf, V_inf, Rho_inf, AR_init, a_inf, Gamma_inf, flag, mu_inf,A2,v2, M2, R,P2,P2_0,h,d2,T2_0, T_h_o, n_modules, T_h_i, T_c_i, T_c_o, T_mean_h, T_mean_c, Q_tot, C_h, C_c, C_star,Pr_air, Pr_coolant, M_dot_coolant, M_dot_2, M_dot_3, T2,P3,P3_0,T3,T3_0,v3,M3,A4,P4,P4_0, T4, T4_0, v4, M4,P5,P5_0,M5,T5,T5_0,v5,P6,P6_0,M6,T6,T6_0,v6, M1, T1, T1_0, Rho_1, v1, P1_0, dia_prop);                 

p_loop = P_inf;
diff_P = P6-p_loop; 
fpr = fpr_init;

 while abs(diff_P) > 10  
    
    [J_val, Cp_val, eff_val, D, J, Cp, V_ind_prop, dp_tot_prop, T_net, dia_prop, M1, P1, P1_0, T1, T1_0, Rho_1, v1,Re_1] = propeller(flight_phase, V_inf, Rho_inf, P_inf, a_inf, T_inf, mu_inf);
    [M3, P3, P3_0, A3, T3, T3_0, Rho_3, v3, M2, P2, P2_0, A2, T2, T2_0, Rho_2, v2, L_diffuser, M_dot_2, M_dot_3, d2, d3, Re_2, Re_3] = Diffuser_og_mine_for_7x7 (M_dot_in, d3, M_dot_FOD, M_dot_comp, M1, P1, P1_0, T1, T1_0, Rho_1, v1, R, theta_max_diff, d2_init, A2, AR_init, a_inf, Gamma_inf, flag, mu_inf) ;                
    [dp_coolant,d_h_air, M_dot_4, b_t_air, b_t_coolant, dp_hx, N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, NTU, R_tot, v_channel_air, v_channel_coolant,d_h_coolant, A_o_coolant, A_o_air, Re_air, Re_coolant, h_air, h_coolant, L_solution, UA_unit, v4, P4_0, M4, T4, T4_0, P4, F_drag_hx, M_hx, A4, drag_HX, Q_pred_sol, T_cool_out, Pr_air, f_air, Nu_air, f_hx, inlet_dp] = HX_design1(e,r,hx_theta, counter, fpr,A3,v3, R,P3,P3_0, h, d3, T_h_o, n_modules, T_h_i, T_c_i, T_c_o, T_mean_h, T_mean_c, Q_tot, C_h, C_c, C_star,Pr_air, Pr_coolant, M_dot_coolant, M_dot_3, T3);
    
    if fan == "ON"
        fpr_init = fan_fpr_simple(M_dot_4, M_dot_des, FPR_des);
    else 
        fpr_init = 1;
    end
    
    [M5,T5_0,T5,P5,P5_0,v5,A5,dp_fan,P_shaft,M_dot_5] = fan_backup (fpr_init,d3,P4,P4_0,T4,T4_0,v4,M4,A4,R,M_dot_4); 
    
    [A6, M6, L_nozzle, P6, P6_0, T6, T6_0, v6, M_dot_6, Re_5,d6] = nozzle_old (AR_noz, flag,d3,P5,P5_0,T5,T5_0,v5,M5,M_dot_4);                 
      
    diff_P = P6-p_loop;
    fprintf("diff p = %f; m_dot_in = %f\n",diff_P, M_dot_in);
    if diff_P < 0
        if diff_P < -3000
            M_dot_in = M_dot_in-(0.1*M_dot_3);
        elseif diff_P > -3000 && diff_P < -1000
            M_dot_in = M_dot_in-(0.01*M_dot_3);
        elseif diff_P > -1000 && diff_P < -200
            M_dot_in = M_dot_in-(0.005*M_dot_3);
        else
            M_dot_in = M_dot_in-(0.0005*M_dot_3);
        end
    elseif diff_P > 0
        if diff_P > 1000
            M_dot_in = M_dot_in+(0.01*M_dot_3);
        elseif diff_P < 1000 && diff_P > 200
            M_dot_in = M_dot_in+(0.005*M_dot_3);
        else
            M_dot_in = M_dot_in+(0.0005*M_dot_3);
        end
    else
        break 
    end
    
    if fan == "OFF"
        if M_dot_in > m_dot_streamtube
            error("Mass flow rate into diffuser is not sufficient; Puller fan is needed;")
        end
    end



    fprintf('P1 = %.2f, P6 = %.2f , diff_P = %.6e, m_dot = %f\n', P1, P6, diff_P,M_dot_in);
 end




fprintf('Converged: |P1 - P6| = %.6e at P_1 = %.2f\n', abs(diff_P), P1);
fprintf('Drag HX: %.2f, %.2f', F_drag_hx, drag_HX)
