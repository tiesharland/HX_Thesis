function drag = Copy_of_tmsObjective(x)
%function P_shaft = Copy_of_tmsObjective(x)   
%%----------------------------------------------------------    INPUTS   -----------------------------------------------------------------------%%

% d2 = x.fr_dia;
d2 = 0.45;
e = x.air_side;
r = x.cool_side;
% r = 4;
AR_noz = 0.36;
AR_diff = 3.8;
% AR_diff = x.AR_diff;
% AR_noz = x.AR_noz;
% N_segments = x.N_segments;
N_segments = 50;

A2               = pi*d2*d2/4;
A_frontal_target = A2*AR_diff;
d3_target        = sqrt(A_frontal_target*4/pi);
AR_diff_arr      = (d3_target/d2)^2;
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

%%%--- MASS FLOW LEAKS ---%%%

M_dot_FOD  = 0;
M_dot_comp = 0;

%%%--- RUN MODEL ---%%%

tol_T = 0; % K
use_DNS = false;

try
    results = run_disc_model_fwdpass(use_DNS, N_segments, e, r, hx_theta, fan, fpr_init, ...
        M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
        d2, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, tol_T);

    drag = results.drag_tot;
    % P_shaft = results.P_shaft;

  disp(class(x))
  whos x
catch ME
    warning('tmsObjective failed');
    fprintf('Error identifier: %s\n', ME.identifier);
    fprintf('Error in: %s, line %d\n', ME.stack(1).name, ME.stack(1).line);
    drag = 1e+03;
    %P_shaft = 1e+05;
    % M_hx = 1e+05;

    return
end


end