function [A6, M6, L_nozzle, P6, P6_0, T6, T6_0, v6, M_dot_6, Re_5,d6,mdot_choked] = nozzle_old (AR_noz, flag,d3,P5,P5_0,T5,T5_0,v5,M5,M_dot_4)                 

  d5 = d3(end);
  d6 = d5*sqrt(AR_noz);
  theta_max_noz = 20; %% GET A GOOD REFERENCE
  Rho_5 = rho_air(P5,T5);
  flag = "Ideal";
  M_dot_5 = M_dot_4;
  %fprintf("Mdot5 = %f\n",M_dot_5);
  L_nozzle = (d5-d6)/(2*tand(theta_max_noz));
  R = 287;

  % loss coefficient estimation based on idelchik pg 375 for alpha 15-40
  % degrees from table 
  AR_k_noz = [0.64, 0.45, 0.39, 0.25, 0.15, 0.1];
  k_arr_noz = [0.040, 0.050, 0.046, 0.044, 0.044, 0.05];
  p_noz = polyfit(AR_k_noz, k_arr_noz, 3);
  k_val_noz = polyval(p_noz, AR_noz);
  fprintf("k_val_noz = %f\n",k_val_noz);
  % k_val_noz = 0;

  %k_val_noz = 0;
  P6_0 = P5_0-((0.5*Rho_5*v5*v5)*k_val_noz);
  T6_0 = T5_0;

  gamma5 = gamma_air(cp_air(T5));
  x = @(M6) (1+(((gamma5-1)*M6^2)/2));
  y = gamma5/(gamma5-1);

  A5 = pi*d5*d5/4;
  A6 = AR_noz*A5; 
  % A_star = M_dot_5/(P5_0*sqrt(gamma_air(cp_air(T5))/(R*T5_0))*(2/(gamma_air(cp_air(T5))+1))^((gamma_air(cp_air(T5))+1)/(2*(gamma_air(cp_air(T5))-1))));
  % fprintf("A6/A_star = %f\n",A6/A_star);

  % Computation of the ideal expansion thru a diffuser i.e no friction effects 
  
   % ---- Choked mass flow ----
    mdot_choked = A6 * P6_0 * sqrt(gamma5/(R*T6_0)) * ...
        (2/(gamma5+1))^((gamma5+1)/(2*(gamma5-1)));
   

    
    if M_dot_5 >= mdot_choked
        M6 = 1;
    else
        f = @(M) M_dot_5 - A6 * P6_0 * sqrt(gamma5/(R*T6_0)) .* M .* (1 + (gamma5-1)/2 .* M.^2).^(-(gamma5+1)/(2*(gamma5-1)));
        M6 = fzero(f, [1e-4, 0.999]);
    end

    T6 = T6_0 / (1 + (gamma5-1)/2 * M6^2);
    gamma6 = gamma_air(cp_air(T6));
    P6 = P6_0 / (1 + (gamma6-1)/2 * M6^2)^(gamma6/(gamma6-1));
    rho6 = P6 / (R*T6);
    v6 = M6 * sqrt(gamma6*R*T6);
    M_dot_6 = rho6*A6*v6;
    Re_5 = Rho_5*A5*v5;
    fprintf("P6 = %f\n",P6);


end
