%Calculating the pressure, density and temperature of the free-stream%
function [p_inf, t_inf, rho_inf, gamma_inf, Cp_inf, mu_inf, P_inf_tot, M_inf, a_inf] = Freestream (h,p11,R,t11,V_inf) 

if (h<11000)
    p_inf =  101325 * (1 - 2.25577e-05 * h)^5.2559; 
    t_inf = (15.04-(0.00649*h))+273.15;
else
    p_inf = p11*exp(-9.8*(h-11000)/(R*t11));
    t_inf = 216.65;
end

rho_inf = p_inf/(R*t_inf);
Cp_inf = 1005 + 0.1*(t_inf-288) - 0.0002*(t_inf-288).^2; %check!!
Cp_inf_check = 1.9327E-10*t_inf^4- 7.9999E-07*t_inf^3+ 1.1407E-03*t_inf^2- 4.4890E-01*t_inf + 1.0575E+03;
gamma_inf = Cp_inf/(Cp_inf-287);
mu_inf = mu_air(t_inf); %4E-29*h^5 - 6E-24*h^4 + 1E-19*h^3 + 1E-14*h^2 - 4E-10*h + 2E-05;  %Digitized from reference plot on oneNote

P_inf_tot = p_inf+(0.5*rho_inf*V_inf*V_inf); %Total free stream pressure - Move to freestream function
a_inf = sqrt(gamma_inf*R*t_inf); %m/s, speed of sound - Move to freestream function
M_inf = V_inf/a_inf; %Freestream Mach number - Move to freestream function



end


