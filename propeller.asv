 function [J_val, Cp_val, eff_val, D, J, Cp, V_ind_prop, dp_tot_prop, T_net, dia_prop, M1, P1, P1_0, T1, T1_0, Rho_1, v1, Re_1] = propeller (flight_phase, V_inf, Rho_inf, P_inf, a_inf, T_inf, mu_inf)

%% READ VALUES FROM EXCEL SHEET %%

D = readtable("C:\Users\prasannamuthuk\Downloads\propellerEfficiency.xlsx",'Sheet','D'); 
J = readtable("C:\Users\prasannamuthuk\Downloads\propellerEfficiency.xlsx",'Sheet','J');
Cp = readtable("C:\Users\prasannamuthuk\Downloads\propellerEfficiency.xlsx",'Sheet','CP');
Ct = readtable("C:\Users\prasannamuthuk\Downloads\propellerEfficiency.xlsx",'Sheet','CT');
eff_1 = readtable("C:\Users\prasannamuthuk\Downloads\propellerEfficiency.xlsx",'Sheet','efficiency_J_CP');
eff_2 = readtable("C:\Users\prasannamuthuk\Downloads\propellerEfficiency.xlsx",'Sheet','efficiency_J_CT');


%% INPUTS %%
n_rpm = 1200; %1200total 
n_rps = n_rpm/60;
dia_prop = 3.96;

%flight_phase = flight_phase;
%fprintf("flight phase = %f\n",flight_phase);

if (flight_phase == 1)
    P_shaft = 1705*1000; %changed from 1775, based on MTO case of PR_HPS_HPS_0045_IssNC - from James
elseif (flight_phase == 3)
    P_shaft = 1475*1000; %PR_HPS_HPS_0044_IssNC - from James
elseif (flight_phase == 2)
    P_shaft = 1475*1000; %at TOC : PR_HPS_HPS_0044_IssNC - from James
end



%% EFFICIENCY VALUE CALCULATION FROM MAP %%

J_val = V_inf/(n_rps*dia_prop);
% fprintf("J_val = %f\n",J_val);
Cp_val = P_shaft/(Rho_inf * n_rps^3 * dia_prop^5);
% fprintf("cp_val = %f\n",Cp_val);

col_val = 0;
row_val = 0;

for i = 1:200
    J_check = round(J{1,i}, 4);
    %fprintf("J_check = %.4f\n", J_check);

    if abs(J_check - J_val) < 1e-2
        col_val = i;
        %fprintf("col = %.4f\n", col_val);
        break; 
    end

end

for m = 1:200
    Cp_check = round(Cp{m,1}, 4);
    %fprintf("cp_check = %.4f\n", Cp_check);

    if abs(Cp_check - Cp_val) < 5e-03
        row_val = m;
        break; 
    end

end

eff_val = eff_1{row_val,col_val};
%% Calculating propeller slipstream characteristics %%

f = @(vi) 2*Rho_inf*(pi/4)*dia_prop*dia_prop*vi*(V_inf+vi).^2 - (P_shaft*eff_val);
V_ind_prop = fzero(f, V_inf);
v1= V_inf+V_ind_prop;
% fprintf("v1 = %f\n",v1);
T_net = (P_shaft*eff_val)/(v1);
fprintf("Induced velocity = %f\n",V_ind_prop);
fprintf("Thrust = %f\n",T_net);

 dp_tot_prop = T_net/(pi*dia_prop*dia_prop/4);
 fprintf("eff_val = %f ;T_net = %f ; dp_tot_prop = %f\n",eff_val,T_net,dp_tot_prop);

 P1_minus = P_inf - (0.5*Rho_inf*((v1^2)-V_inf^2));
 P1 = P1_minus+dp_tot_prop;
 fprintf("P1 = %f\n",P1);
 M1 = v1/a_inf;
 fprintf("M1 = %f\n",M1);
 fprintf("v1 = %f\n",v1);
 P1_0 = P1*((1+((gamma_air(cp_air(T_inf))-1).*M1.^2/2)).^((gamma_air(cp_air(T_inf)))./(gamma_air(cp_air(T_inf))-1)));
 Rho_1 = rho_air(P1,T_inf);
 T1 = T_inf;
 T1_0 = T_inf*(1+((gamma_air(cp_air(T_inf))-1)*M1.^2)/2);
 Re_1 = Rho_1*v1*dia_prop/mu_air(T1);
 

end
