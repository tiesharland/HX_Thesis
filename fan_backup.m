function [M5,T5_0,T5,P5,P5_0,v5,A5,dp_fan,P_shaft,M_dot_5] = fan_backup (fpr, d3,P4,P4_0,T4,T4_0,v4,M4,A4,R,M_dot_4)
    R = 287;

    fan_eff = 0.9; 
    A5 = A4;
    P5_0 = fpr*P4_0;

    %Estimation of total temperature increase%
    T5_0_ideal = T4_0 * (fpr)^((gamma_air(cp_air(T4))-1)/gamma_air(cp_air(T4)));
    T5_0 = T4_0 + ((T5_0_ideal-T4_0)/fan_eff);

    M5 = M4;
    %fprintf("M5 = %f\n",M5);

    T5 = T5_0/(1+((gamma_air(cp_air(T4))-1)*M5.^2)/2);
    P5 = P5_0/((1+((gamma_air(cp_air(T5))-1).*M5.^2/2)).^((gamma_air(cp_air(T5)))./(gamma_air(cp_air(T5))-1)));
    v5 = M5*sqrt(gamma_air(cp_air(T5))*R*T5);
    fprintf("P5 = %.2f\n",P5);

    m_nozz = A5*P5_0*sqrt(gamma_air(cp_air(T5))/(R*T5_0))*(2/(gamma_air(cp_air(T5))+1))^((gamma_air(cp_air(T5))+1)/(2*(gamma_air(cp_air(T5))-1)));

    dp_fan = P5-P4;
    %fprintf("dp_fan_term = %.2f Pa \n", dp_fan)
    %M_dot_5 = A5 * P5_0 * sqrt(gamma_air(cp_air(T5))/(R*T5_0))* M5 * (1 + (gamma_air(cp_air(T5))-1)/2 * M5^2)^ (-(gamma_air(cp_air(T5))+1)/(2*(gamma_air(cp_air(T5))-1)));
    %M_dot_5 = M_dot_4;
    M_dot_5 = min(M_dot_4, m_nozz);
    fprintf("mdot5 = %f\n",M_dot_5)
    P_fan = dp_fan*M_dot_5/rho_air(P5,T5);
    fprintf("mdot5 = %f\n",M_dot_5)
    P_shaft = P_fan/fan_eff;
    fprintf("Fan power = %f kW\n\n\n",P_shaft/1000);


%end

