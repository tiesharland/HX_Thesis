clear all
clc

counter = 1;

e = 6;
r = 4.2;
d2_init = 0.52; %0.45808;
AR_diff = 3.8;
AR_noz  = 0.36;
A2      = pi*d2_init*d2_init/4;
A_frontal_target = A2*AR_diff;
d3_target        = sqrt(A_frontal_target*4/pi);
AR_diff_arr      = (d3_target/d2_init)^2;
n_modules        = 2;
M_dot_coolant    = 44.4;

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
fpr_init = 1;

fan  = "OFF";
disc = "ON";
N_segments = 50;

%%%--- MASS FLOW LEAKS ---%%%

M_dot_FOD  = 0;
M_dot_comp = 0;

%%%--- DIFFUSER INPUTS ---%%%

AR_init = AR_diff_arr;

%%%--- RUN MODEL ---%%%

tol_T = 0; % K
use_DNS = true;

results_D = run_disc_model_fwdpass(use_DNS, N_segments, e, r, hx_theta, fan, fpr_init, ...
    M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
    d2_init, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, tol_T);


%%%--- SUMMARY PRINTS ---%%%

fprintf("\n=== HX RESULTS ===\n");
fprintf("HX length          = %.4f m\n",   results_D.L_solution);
fprintf("HX pressure drop   = %.2f Pa\n",  results_D.dp_hx);
fprintf("Air outlet temp    = %.2f K\n",   results_D.T4);
fprintf("HX mass            = %.2f kg\n",  results_D.M_hx);
fprintf("HX drag         = %.2f N\n",   results_D.drag_HX);
fprintf("P6 - P_inf         = %.2f Pa\n",  results_D.P6 - results_D.P_inf);
fprintf("Coolant outlet temp = %.2f K\n",   results_D.T_cool_out);
fprintf("Coolant pressure drop (avg) = %.2f Pa\n",   results_D.dp_coolant);
fprintf("Cooling power  = %.3f MW\n",   results_D.Q_pred_solution/1e6);
fprintf("Nozzle exit P6     = %.2f Pa\n",  results_D.P6);
fprintf("Ambient P_inf      = %.2f Pa\n",  results_D.P_inf);
fprintf("d_h / L            = %f  \n", results_D.d_h_air/results_D.L_solution);

%%%--- PLOTS ---%%%

plot_stations(results_D);

plot_HX(results_D);