function [M3, P3, P3_0, A3, T3, T3_0, Rho_3, v3, M2, P2, P2_0, A2, T2, T2_0, Rho_2, v2, L_diffuser, M_dot_in, M_dot_out, d2, d3, Re_2, Re_3] = Diffuser_og_mine_for_7x7 (M_dot_in, d3, M_dot_FOD, M_dot_comp, M1, P1, P1_0, T1, T1_0, Rho_1, v1, R, theta_max_diff, d2, A2, AR_init, a_inf, Gamma_inf, flag, mu_inf)                 
    
  AR = AR_init;

  % loss coefficient estimation based on idelchik figure 5.12 from edition at 12 degrees
  AR_k = [2, 3, 4, 6];
  k_arr = [0.28 0.21 0.19 0.185]; %at 12 degrees
  %k_arr = [0.32 0.26 0.26 0.26]; % at 15 degrees
  %k_arr = [0.39 0.39 0.39 0.39]; %at 20 degrees
  p = polyfit(AR_k, k_arr, 2);
  k_val = polyval(p, AR);
  fprintf("k_val = %f\n",k_val);
  %k_val = 0;


  % Computation of the ideal expansion thru a diffuser i.e no friction effects % 
  T2_0 = T1_0;
  T2 = T1;

  P2_0 = P1_0*0.99; %assuming an inlet loss of 1 percent 
  %fprintf("M1 = %f\n",M1)
  

  v2 = M_dot_in/(Rho_1*A2);
  %fprintf("v1 , v2 = %f %f\n",v1,v2)
  gamma2 = gamma_air(cp_air(T2));
  x = @(M2) (1+(((gamma2-1)*M2^2)/2));
  y = gamma2/(gamma2-1);

  f_M2 = @(M2) ((P2_0)/(R*T2_0))*x(M2).^(-(y-1))*A2*M2*sqrt(gamma2*R*T2_0/x(M2))-M_dot_in;
  M2_low = 0.01;
  M2_high = 0.99;

  M2 = fzero(f_M2, [M2_low, M2_high]);
  %fprintf("M2 = %.2f\n",M2);
 
  P2=P2_0/(((1+((gamma_air(cp_air(T2))-1).*M2.^2/2)).^((gamma_air(cp_air(T2)))./(gamma_air(cp_air(T2))-1))));

  Rho_2 = P2 / (R * T2);
  a_2 = sqrt(gamma_air(cp_air(T1))*287*T1);
  v2 = a_2*M2;
  fprintf("P2 = %f\n",P2)

  A3 = pi*d3*d3/4;
  L_diffuser = (d3-d2)/(2*tand(theta_max_diff));
  %fprintf("L-diff = %f \n",L_diffuser);
  Re_2 = Rho_2*v2*d2/mu_air(T2); 
  A2_Athroat = (1/M2) * ((2/(Gamma_inf+1)) * (1 + (Gamma_inf-1)/2 * M2^2))^((Gamma_inf+1)/(2*(Gamma_inf-1)));
  A_throat = A2/A2_Athroat;
  
  M3_ideal = fzero(@(M3) (1/M3) * ((2/(Gamma_inf+1)) * (1 + (Gamma_inf-1)/2 * M3^2))^((Gamma_inf+1)/(2*(Gamma_inf-1))) - (A3/A_throat), 0.2);  
  
  if (flag=="Ideal")
      M3=M3_ideal;
  else
      % Computation of expansion thru a diffuser with friction %
      f= 0.005; %%Scope for improvement - this value holds good for Re > 10^5, surface roughness - 0.001*D; 
      d_x = @(x) d2+(d2*(sqrt(AR)-1)*(x/L_diffuser));
      A = @(x) A2 * (1 + ((AR-1)*x/L_diffuser));      
      dAdx = @(x) ((AR-1)/L_diffuser)*A2; 
      dMdx = @(x,M) (M ./ (1 - M.^2)) .* ( -( (1 + 0.5*(Gamma_inf-1).*M.^2) ./ A(x) ) .* dAdx(x) + ( (0.5*Gamma_inf.*M.^2.*f) ./ d_x(x) ) .* (2 + (Gamma_inf-1).*M.^2) ./ 2 );

      tspan = linspace(0,L_diffuser,100);
      options = odeset('Events', @(x,M) choke_event(x,M));
      [~, M_sol] = ode45(dMdx,tspan, M2,options); 
      M3=M_sol(end);
  end

  fprintf("M3 = %.2f\n",M3);

T3_0 = T2_0;
T3 = T3_0/(1+((gamma_air(cp_air(T2))-1)/2)*M3^2);
P3_0 = P2_0-((0.5*Rho_2*v2*v2)*k_val); 
P3 = P3_0/((1+((gamma_air(cp_air(T3))-1).*M3.^2/2)).^((gamma_air(cp_air(T3)))./(gamma_air(cp_air(T3))-1)));
Rho_3 = P3/(R*T3); 
a_3 = sqrt(gamma_air(cp_air(T3))*287*T3);
v3 = a_3*M3;
fprintf("P3 = %.2f\n",P3);
M_dot_out = (Rho_3*A3*v3(end))-(M_dot_in*M_dot_FOD)-(M_dot_in*M_dot_comp);
Re_3 = Rho_3*v3*d3/mu_air(T3);


end






  