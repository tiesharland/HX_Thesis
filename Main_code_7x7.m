clear all
clc

% w = warning('query','last');
% id = w.identifier;
% warning('off',id);
drag_tot = zeros(5,5);
l_hx_arr = zeros(5,5);
M_hx_arr = zeros(5,5);

f = 0; 

counter=1;


for i = 1:1

    for j = 1:1

        d2_init = [0.46]; %m, diffuser inlet initial guess
        A2 = pi*d2_init*d2_init/4;
        A_frontal_target = A2*4; %this number multiplied to A2 is your diffuser AR
        d3_target = sqrt(A_frontal_target.*4/pi);
        AR_diff_arr = (d3_target./d2_init).^2;
        n_modules = 2; %ducts per nacelle

%%%--------------------------------INPUTS - FREESTREAM CALCS------------------------------------%%%

flight_phase = 3; % 1 - take-off; 2 - top of climb; 3 - cruise;

if (flight_phase == 1)
    h = 5; %m
    T_in_fc = 85+273; %K; HX outlet coolant temperature from James'email
    T_out_fc = 105+273; %K; HX inlet coolant temperature from James'email
    Q_tot = 2.38*1e+06; %W, total heat produced/cooling power reqd
    V_inf = 2;
    k = 5;
elseif (flight_phase == 3)
    h = 4876; %m
    T_in_fc = 70+273; %K; HX outlet coolant temperature from James'email
    T_out_fc = 85+273; %K; HX intlet coolant temperature from James'email
    V_inf = 128;
    Q_tot = 2.25e+06; %W, total heat produced/cooling power reqd
    k=1;
elseif (flight_phase == 2)
    h = 4876; %m
    T_in_fc = 70+273; %K; HX outlet coolant temperature from James'email
    T_out_fc = 85+273; %K; HX inlet coolant temperature from James'email
    Q_tot = 2.25*1e+06; %W, total heat produced/cooling power reqd
    V_inf = 128;
    k=1;

end

ref_pressure = 101325; % Pa
Lapse_rate = 0.0065; % K/m, lapse rate
p11 = 22632; % Pa, constant value for calculation of pressure > 11km altitude
t11 = 216.65; % K, constant value for calculation of pressure > 11km altitude
R = 287; % J/kg.K, gas constant
gamma = 1.4;
hx_theta = 60; %hx inclination angle
counter = 1; 


fan = "OFF";
disc = "ON";

%% ----------------------------------MASS FLOW LEAKS---------------------------------------%%

M_dot_FOD = 0; % Assumption - 5% of the total ingested mass flow goes to the foreign object Debris 
M_dot_comp = 0; % Assumption - 5% of the total ingested mass flow goes to the compressor for the fuel cell

%%%--------------------------------INPUTS - DIFFUSER--------------------------------------------%%%

A2_init = pi*(d2_init^2)/4; %m^2
AR_init = AR_diff_arr; %initial assumption for area ratio (A2/A1)

%% ----------------------------------------------Freestream calculations -------------------------------------------- %%

[P_inf, T_inf, Rho_inf, Gamma_inf, Cp_inf, mu_inf, P_inf_tot, M_inf, a_inf] = Freestream (h,p11,R,t11,V_inf);

P_inf_tot = P_inf+(0.5*Rho_inf*V_inf*V_inf); %Total free stream pressure - Move to freestream function
a_inf = sqrt(Gamma_inf*R*T_inf); %m/s, speed of sound - Move to freestream function
M_inf = V_inf/a_inf; %Freestream Mach number - Move to freestream function


%% ----------------------------------------------- Propeller ----------------------------------------------%%%

[J_val, Cp_val, eff_val, D, J, Cp, V_ind_prop, dp_tot_prop, T_net, dia_prop, M1, P1, P1_0, T1, T1_0, Rho_1, v1,Re_1] = propeller(flight_phase, V_inf, Rho_inf, P_inf, a_inf, T_inf, mu_inf);


%% ----------------------------------------------- Diffuser ----------------------------------------------%%%
theta_max_diff = 15; %diffuser expansion half angle
flag = "Friction"; %"Ideal" or "Friction"
d3 = d3_target;
disp(d3);
m_dot_streamtube = (Rho_1*v1*A2);

if fan == "OFF"
    M_dot_in = m_dot_streamtube;
else
    M_dot_in = m_dot_streamtube*k;
end

[M3, P3, P3_0, A3, T3, T3_0, Rho_3, v3, M2, P2, P2_0, A2, T2, T2_0, Rho_2, v2, L_diffuser, M_dot_2, M_dot_3, d2, d3, Re_2, Re_3] = Diffuser_og_mine_for_7x7 (M_dot_in, d3, M_dot_FOD, M_dot_comp, M1, P1, P1_0, T1, T1_0, Rho_1, v1, R, theta_max_diff, d2_init, A2, AR_init, a_inf, Gamma_inf, flag, mu_inf) ;                



%% --------------------------------INPUTS - Heat exchanger general -------------------------------------- %%


delta_T = T_out_fc-T_in_fc;
M_dot_coolant_max = 44.4; %UPDATED VALUES FROM SAM on 29/01/2026; kg/s; from excel sheet "02. Air parameters TUD June 2025 to share - issue A" 
M_dot_coolant_min = 29; %kg/s; from excel sheet "02. Air parameters TUD June 2025 to share - issue A" 
M_dot_coolant = M_dot_coolant_max;

%%%---------------------------------------- Mean temperature calculations -----------------------------------------------%%%

[T_mean_h, T_mean_c, T_c_i, T_c_o, T_h_i, T_h_o, C_h, C_c, C_star] = HX_deltaT (T3, M_dot_3, M_dot_coolant,T_out_fc,T_in_fc);

%%%--------------------------------------- Prandtl number calculations --------------------------------------------------%%%

Pr_coolant = mu_EG(T_mean_h)*cp_EG_50_50(T_mean_h)/k_EG(T_mean_h);
Pr_air = mu_air(T_mean_c)*cp_air(T_mean_c)/k_air(T_mean_c);

%%%----------------------------------------------- Sizing function ----------------------------------------------%%%
x=1;
e = 8; %this is a multiplier fed into HX sizing function to indicate air-side fin size (1 = 1mm)
r = 4.2; %this is a multiplier fed into HX sizing function to indicate coolant-side fin size (1 = 1mm)
[dp_coolant,d_h_air, M_dot_4, b_t_air, b_t_coolant, dp_hx, N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, NTU, R_tot, v_channel_air, v_channel_coolant,d_h_coolant, A_o_coolant, A_o_air, Re_air, Re_coolant, h_air, h_coolant, L_solution, UA_unit, v4, P4_0, M4, T4, T4_0, P4, F_drag_hx, M_hx, A4, drag_HX, Q_pred_sol, T_cool_out, Pr_air, f_air, Nu_air, inlet_dp] = HX_design1(e,r,hx_theta, counter, x ,A3,v3, R,P3,P3_0, h, d3, T_h_o, n_modules, T_h_i, T_c_i, T_c_o, T_mean_h, T_mean_c, Q_tot, C_h, C_c, C_star,Pr_air, Pr_coolant, M_dot_coolant, M_dot_3, T3);         
coolant_pres_drop(i,j) = dp_coolant;
HX_dp (i,j) = dp_hx;

%% ----------------------------------------------- Puller fan ---------------------------------------------- %%
M_dot_des = rho_air(P4,T4)*A4*M4*sqrt(gamma_air(cp_air(T4))*R*T4); %design point mass flow rate
FPR_des = 1.01; %design point FPR
if fan == "ON"
    fpr_init = fan_fpr_simple(M_dot_4, M_dot_des, FPR_des);
    fprintf("Fan pressure ratio = %f\n",fpr_init);
else 
    fpr_init = 1;
end

[M5,T5_0,T5,P5,P5_0,v5,A5,dp_fan,P_shaft,M_dot_5] = fan_backup (fpr_init,d3,P4,P4_0,T4,T4_0,v4,M4,A4,R,M_dot_4); 

%% ----------------------------------------------- Nozzle design ----------------------------------------------%%%

AR_noz = 0.33; % Nozzle area ratio
[A6, M6, L_nozzle, P6, P6_0, T6, T6_0, v6, M_dot_6, Re_5,d6,mdot_choked] = nozzle_old (AR_noz, flag,d3,P5,P5_0,T5,T5_0,v5,M5,M_dot_4);           


%% ----------------------------------------------- Pressure difference correction ----------------------------------------------%%

% Without HX_deltaT in iteration
[P1,P2,P3,P4,P5,P6,P1_0,P2_0,P3_0,P4_0,P5_0,P6_0, M1, M2, M3, M4, M5, M6,  v1,v2,v3,v4,v5,v6, T1, T1_0, T2, T2_0, T3,T3_0,T4, T4_0, T5, T5_0, T6, T6_0, M_dot_2, M_dot_3, M_dot_4, M_dot_5, M_dot_6,P_shaft, M_hx, dp_coolant, L_solution, Q_pred_sol, T_cool_out, dp_hx, drag_HX, h_air, v_channel_air, Re_air, Pr_air, f_air, Nu_air, f_hx, inlet_dp] = pressure_eq_mDot_adjust_9thJan_noExhaust(e,r,m_dot_streamtube, fan, M_dot_des, FPR_des, hx_theta,fpr_init, M_dot_in, flight_phase, theta_max_diff,d3,M_dot_FOD, M_dot_comp,counter, A2_init, d2_init,AR_noz, dp_tot_prop, P1,P_inf, P_inf_tot, T_inf, V_inf, Rho_inf, AR_init, a_inf, Gamma_inf, flag, mu_inf,A2,v2, M2, R,P2,P2_0,h,d2,T2_0, T_h_o, n_modules, T_h_i, T_c_i, T_c_o, T_mean_h, T_mean_c, Q_tot, C_h, C_c, C_star,Pr_air, Pr_coolant, M_dot_coolant, M_dot_2, M_dot_3, T2,P3,P3_0,T3,T3_0,v3,M3,A4,P4,P4_0, T4, T4_0, v4, M4,P5,P5_0,M5,T5,T5_0,v5,P6,P6_0,M6,T6,T6_0,v6, M1, T1, T1_0, Rho_1, v1, P1_0, dia_prop);
% % With HX_deltaT in iteration
% [P1,P2,P3,P4,P5,P6,P1_0,P2_0,P3_0,P4_0,P5_0,P6_0, M1, M2, M3, M4, M5, M6,  v1,v2,v3,v4,v5,v6, T1, T1_0, T2, T2_0, T3,T3_0,T4, T4_0, T5, T5_0, T6, T6_0, M_dot_2, M_dot_3, M_dot_4, M_dot_6,P_shaft, M_hx, dp_coolant] = pressure_iteration_deltaT(e,r,m_dot_streamtube, fan, M_dot_des, FPR_des, hx_theta,fpr_init, M_dot_in, flight_phase, theta_max_diff,d3,M_dot_FOD, M_dot_comp,counter, A2_init, d2_init,AR_noz, dp_tot_prop, P1,P_inf, P_inf_tot, T_inf, V_inf, Rho_inf, AR_init, a_inf, Gamma_inf, flag, mu_inf,A2,v2, M2, R,P2,P2_0,h,d2,T2_0, n_modules, Q_tot, Pr_air, Pr_coolant, T_out_fc,T_in_fc, M_dot_coolant, M_dot_2, M_dot_3, T2,P3,P3_0,T3,T3_0,v3,M3,A4,P4,P4_0, T4, T4_0, v4, M4,P5,P5_0,M5,T5,T5_0,v5,P6,P6_0,M6,T6,T6_0,v6, M1, T1, T1_0, Rho_1, v1, P1_0, dia_prop);                 

%% ----------------------------------------------- State variables - Plot ----------------------------------------------%%

plot1 (AR_init, AR_noz, fpr_init, P_inf, P1, P2, P3, P4, P5,P6, T_inf,T1, T2, T3, T4, T5, T6, P_inf_tot, P1_0, P2_0, P3_0, P4_0, P5_0, P6_0, V_inf, v1, v2, v3, v4, v5, v6, v_channel_air, M_inf, M1, M2, M3, M4, M5, M6);


if fpr_init == 1
    fprintf("The puller fan is not present\n")
else 
    fprintf("The puller fan is present\n")
end


m_dot_spill = m_dot_streamtube-M_dot_2;
fprintf("Spillage mass flow rate = %f\n",m_dot_spill);
drag_tot(i,j) = (m_dot_streamtube*v1) - (M_dot_6*v6) - (M_dot_2*M_dot_FOD*v3) + ((P2-P_inf)*A2) - ((P6-P_inf)*A6);




fprintf("Total drag = %f\n",drag_tot);

    end
end

