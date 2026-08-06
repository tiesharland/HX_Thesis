function [dp_coolant_loop,d_h_air, M_dot_4, b_t_air, b_t_coolant, dp_hx, N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, NTU, R_tot, v_channel_air, v_channel_coolant,d_h_coolant, A_o_coolant, A_o_air, Re_air, Re_coolant, h_air, h_coolant, L_solution, UA_unit, v4, P4_0, M4, T4, T4_0, P4, F_drag_hx, M_hx, A4, drag_HX, Q_pred_solution, T_cool_out, Pr_air, f_air, Nu_air, f_hx, inlet_dp] = HX_design1(e,r,hx_theta, counter, ~,A3,v3, R,P3,P3_0, h, d3, T_h_o, n_modules, T_h_i, T_c_i, T_c_o, T_mean_h, T_mean_c, Q_tot, C_h, C_c, C_star,Pr_air, Pr_coolant, M_dot_coolant, M_dot_3, T3)

%% Inputs %%

A_frontal = pi*d3*d3/(4*cosd(hx_theta)); %m^2; Total frontal area % 0.3969; %
fprintf("HX frontal area = %f\n",A_frontal);
fprintf("e = %f\n",e);

m_dot_air = M_dot_3;
%fprintf("m_dot_air in hx = %f\n",m_dot_air)
Q_limit = m_dot_air*cp_air(T3)*(T_h_i-T_c_i);
fprintf("Qlim = %f\n",Q_limit);
if Q_limit < Q_tot/n_modules
   error ("HX mass flow rate not sufficient, increase diffuser inlet area or need a puller fan;");
end
m_dot_coolant = M_dot_coolant; %kg/s; Coolant mass flow rate

t_plate = 0.5e-03; %m; thickness of the plate in between consecutive air and coolant passages
t_fin = 0.1e-03; %m; thickness of the fin

b_t_air = 1e-03*e; %1.5e-03; %m; side of the equilateral triangle on air side
b_t_coolant = 1e-03*r; %2e-03; %m; side of the equilateral triangle on coolant side

pitch_fin_air = (0.5*b_t_air)+t_fin; %distance between start of fin and mid point of next fin
pitch_fin_coolant = (0.5*b_t_coolant)+t_fin; %distance between start of fin and mid point of next fin

ht_air = (sqrt(3)/2)*b_t_air; %m; calculation of height of triangular fin
ht_coolant = (sqrt(3)/2)*b_t_coolant;
%fprintf("Height of fin air side; coolant side = %f %f \n",ht_air, ht_coolant);

h_passage_air = ht_air+(t_fin*sqrt(2)); %m; space in between two consecutive plates
h_passage_coolant = ht_coolant+(t_fin*sqrt(2));

%% Calculating N passes %%

b_hx = sqrt(A_frontal); %Assuming the hx frontal area is equal to area of circle at diffuser outlet
h_hx = b_hx;
fprintf("HX frontal area = %f\n",b_hx*h_hx);
N_pass_tot = h_hx/((t_plate+(ht_air+(t_fin/sqrt(2))))+((t_plate+(ht_coolant+(t_fin/sqrt(2)))))); %total number of passages (air+coolant); previous one was without the second term in the denominator - changed: 03/03/26
fprintf("Passages = %f\n",N_pass_tot);
% N_pass_tot_prev = h_hx/(t_plate+(ht_air+(t_fin/sqrt(2)))); %total number of passages (air+coolant); previous one was without the second term in the denominator - changed: 03/03/26
% fprintf("Passages prev = %f\n",N_pass_tot_prev);
N_fin_air = 2*floor(b_hx/((2*t_fin/sqrt(2))+b_t_air)); %number of fins calculated thru total breadth of hx
N_air_pass = floor(N_pass_tot)+1; %number of air passages floor(N_pass_tot*0.5)+1;
N_coolant_pass = N_air_pass-1; %number of coolant passages

%% Hydraulic diameter %%

d_h_air = 4*(0.5*b_t_air*ht_air)/(3*b_t_air);
d_h_coolant = 4*(0.5*b_t_coolant*ht_coolant)/(3*b_t_coolant);

%% Free flow area calculation %%

A_o_air = (0.5*(b_t_air*0.5)*ht_air)+(0.5*((b_t_air*0.5)+(t_fin/sqrt(2)))*(ht_air+(t_fin/sqrt(2))));
A_o_coolant = (0.5*(b_t_coolant*0.5)*ht_coolant)+(0.5*((b_t_coolant*0.5)+(t_fin/sqrt(2)))*(ht_coolant+(t_fin/sqrt(2))));  
%A_o_air = 0.5*(b_t_air - t_fin*sqrt(3))*(ht_air - t_fin);
%A_o_coolant = 0.5*(b_t_coolant - t_fin*sqrt(3))*(ht_coolant - t_fin);

%% Reynolds calculations - air side %%

m_dot_pass_air = m_dot_air/N_air_pass;
%fprintf("md_air = %f\n",m_dot_pass_air);
m_dot_pass_coolant = m_dot_coolant/N_coolant_pass;

v_channel_air = (m_dot_pass_air/(N_fin_air))/(rho_air(P3,T_mean_c)*(A_o_air));

Re_air = rho_air(P3,T_mean_c)*v_channel_air*d_h_air/mu_air(T_mean_c);
%fprintf("Re = %f\n", Re_air);


%% j and f calculations - air side %%

f_air_og = 2.69*(Re_air^(-0.918))*((ht_air/(b_t_air+t_fin))^0.355)*((t_fin/(b_t_air+t_fin))^-0.175);
%fprintf("Old friction factor. for air side = %f\n",f_air_og);
j_air = 0.789*(Re_air^(-1.1218))*((ht_air/(b_t_air+t_fin))^1.235)*((t_fin/(b_t_air+t_fin))^-0.764);
G_air = rho_air(P3,T_mean_h)*v_channel_air;
Nu_air_og = 0.0214*((Re_air^0.8)-100)*(Pr_air^0.4)*(1+(d_h_air/h_hx)^(2/3))*(T_mean_c/T_mean_h)^0.45; %Last term is T_mean_hx/T_wall; T_mean_hx is calculated thru Q/mdotCp for air side; T_wall is assumed to be coolant inlet temp;
%fprintf("Old Nusselt no. for air side = %f\n",Nu_air_og);
%h_air = Nu_air*k_air(T_mean_c)/d_h_air;

%% Effectiveness calculations %% 

    % function epsilon = eps_crossflow_unmixed(NTU, C_star) %effectiveness of fin calculations using NTU
    %     term = (NTU^0.22)*(exp(-C_star*NTU^0.78)-1)/C_star; 
    %     epsilon = (1 - exp(term))*1.5;
    % end

    % function eps_values = solve_eps(NTU, C_star)
    % % Coefficients of cubic: a*x^3 + b*x^2 + c*x + d = 0
    % a = 23.3951;
    % b = -31.0889 + 18.2731*C_star;
    % c = 13.3405 - 14.6761*C_star + 7.25239*(C_star^2);
    % d = -1.16111 + 2.10253*C_star - 1.02985*(C_star^2) ...
    %     - 0.191992*(C_star^3) - NTU;
    % 
    % % Solve cubic
    % roots_all = roots([a b c d]);
    % 
    % % Filter for real roots in [0,1]
    % eps_values = roots_all(imag(roots_all)==0 & real(roots_all)>=0 & real(roots_all)<=1);
    % eps_values = real(eps_values); % ensure numeric type
    % 
    % end

function eps = epsilon_counterflow(NTU, C_star)

    if C_star < 1e-3     % avoid stiffness
        eps = 1 - exp(-NTU); 
    else
        eps = (1 - exp(-NTU .* (1 - C_star))) ./ ...
          (1 - C_star .* exp(-NTU .* (1 - C_star)));
    end
end


%% Heat transfer prediction based on L_hx %%

A_ht_air = 0;
k_fin = 0;
m_air = 0;
eta_fin_air = 0;

  function [Qpred] = Q_pred_L(L) 
    
    %----------------- Air side ------------------%  

    A_p_air = (b_t_air*0.5)+((b_t_air*0.5)+(t_fin/sqrt(2))); %m; air side primary heat transfer surface
    A_f_air = sqrt(ht_air^2 + (b_t_air*0.5)^2)+sqrt(((b_t_air*0.5)+(t_fin/sqrt(2)))^2 + (ht_air+(t_fin/sqrt(2)))^2); %m; air side fin surface area
    %fprintf("Total ht sa = %f\n",A_ht_air)

    A_ht_air_1 = (((A_f_air*N_fin_air*N_air_pass) + (A_p_air*N_fin_air*N_air_pass))*L); %m^2; total heat transfer area
    %fprintf("Total ht sa 1 = %f\n",A_ht_air_1);

    A_ht_air = 3*b_t_air*2*N_fin_air*N_air_pass*L;
    %fprintf("Total ht sa = %f\n",A_ht_air);

    f_air = 0.25*((1.8*log10(Re_air))-1.5)^(-2); %Konakov et al as cited in Mortean et.al 2019
    %fprintf("ff for air side = %f\n",f_air);
    K = T_mean_h/T_mean_c;
    if Re_air < 2300
        Nu_air = 3.66;  % Graetz solution for constant wall temperature
    elseif Re_air < 3000
        w = (Re_air - 2300)/(3000 - 2300);
        Nu_lam = 3.66;
        Nu_turb = (f_air/2)*(Re_air-1000)*Pr_air*(1+(d_h_air/L)^(2/3))*K/(1+(12.7*sqrt(f_air/2)*(Pr_air^(2/3)-1))); %Gnielinski et.al as cited in Mortean et.al 2019
        Nu_air = (1-w)*Nu_lam + w*Nu_turb;
    else
    % Turbulent Gnielinski
        Nu_air = (f_air/2)*(Re_air-1000)*Pr_air*(1+(d_h_air/L)^(2/3))*K/(1+(12.7*sqrt(f_air/2)*(Pr_air^(2/3)-1))); %Gnielinski et.al as cited in Mortean et.al 2019
    end

    %fprintf("New Nusselt no. for air side = %f\n",Nu_air);
    h_air = Nu_air*k_air(T_mean_c)/d_h_air;
    %fprintf("h air = %f\n",h_air);

    k_fin = 237; %W/mK; thermal conductivity of aluminium
    m_air = sqrt(2*h_air/(k_fin*t_fin)); %constant used in fin efficiency calcs.
    eta_fin_air = tanh(m_air*ht_air)/(m_air*ht_air); %fin efficiency

    R_cond = (t_plate/(k_fin*t_plate*(2*N_pass_tot + 2)*b_hx)); %K/W; conductive heat transfer resistance

    R_air = 1 / (h_air * A_ht_air *eta_fin_air); %K.m/W; convective heat transfer resistance air side

    %----------------- Coolant side ------------------%  

    N_fin_coolant = 2*L/(2*t_fin/sqrt(2)+b_t_coolant);

    A_p_coolant = (b_t_coolant*0.5)+((b_t_coolant*0.5)+(t_fin/sqrt(2))); %m^2; coolant side primary heat transfer surface
    A_f_coolant = sqrt(ht_coolant^2 + (b_t_coolant*0.5)^2)+sqrt(((b_t_coolant*0.5)+(t_fin/sqrt(2)))^2 + (ht_coolant+(t_fin/sqrt(2)))^2);  %m^2; coolant side fin surface area
    
    A_ht_coolant = (A_f_coolant*N_fin_coolant*N_coolant_pass) + (A_p_coolant*N_fin_coolant*N_coolant_pass); %m^2; total heat transfer area
    
    v_channel_coolant = (m_dot_pass_coolant/(2*N_fin_coolant))/(rho_EG(T_mean_h)*(A_o_coolant));
    Re_coolant = rho_EG(T_mean_h)*v_channel_coolant*d_h_coolant/mu_EG(T_mean_h);


    %fprintf("Coolant side friction factor = %f  \n",f_coolant);
    

    if Re_coolant < 2300
        % ---------------- LAMINAR FLOW ----------------
   
        Nu_coolant = 4.36;  % constant wall temperature (equilateral triangular passage)
        
        h_coolant = Nu_coolant * k_EG(T_mean_h) / d_h_coolant;
    
    else
        % ---------------- TURBULENT FLOW ----------------
    
        Nu_coolant = 0.023 * Re_coolant^0.8 * Pr_coolant^0.4;
   
        j_coolant = Nu_coolant / (Re_coolant * Pr_coolant^(1/3));
        Cp_c = cp_EG_50_50(T_mean_h);
        G_coolant = rho_EG(T_mean_h)*v_channel_coolant;
        h_coolant = j_coolant * G_coolant * Cp_c / Pr_coolant^(2/3);

    end
    
    m_coolant = sqrt(2*h_coolant/(k_fin*t_fin)); %constant used in fin efficiency calcs.
    eta_fin_coolant = tanh(m_coolant*ht_coolant)/(m_coolant*ht_coolant);
    R_cool = 1 / (h_coolant * A_ht_coolant *eta_fin_coolant); %K/W; convective heat transfer resistance coolant side
    R_tot = R_air + R_cond + R_cool; %K/W; total heat transfer resistance

    UA_unit = 1 / R_tot;

    NTU = (UA_unit)/min(C_c, C_h);
    %fprintf("NTU = %f\n", NTU);
    %x = @(eps) (-1.16111 +(13.3405*eps)+(2.10253*C_star)-(31.0889*eps^2)-(14.6761*eps*C_star)-(1.02985*C_star^2)+(23.3951*eps^3) + (18.2731*eps^2*C_star)+(7.25239*eps*C_star^2)-(0.191992*C_star^3)); %Digitized plot and fitted fro Shah and Sekulic
    %epsilon = eps_crossflow_unmixed(NTU, C_star);   % pick according to mixing case

    epsilon = epsilon_counterflow(NTU, C_star);

    %fprintf("epsilon = %f\n", epsilon);
    if isempty(epsilon)

     % Warning for debugging (remove or comment later)
        warning('HX_design1:solve_eps returned empty for NTU=%g, C_star=%g. Using fallback epsilon.', NTU, C_star);

        % Reasonable fallback approximations (choose one):
        % 1) Simple single-pass approximation (monotonic in NTU):
        epsilon_fallback = 1 - exp(-NTU);           % simple, safe fallback
        % 2) Alternative fallback for small C_star (unbalanced): uncomment if desired
        % epsilon_fallback = 1 - exp(-NTU*(1 - C_star));
        epsilon = min(1, max(0, epsilon_fallback)); % clamp 0..1

    end


    if ~isscalar(epsilon)
         epsilon = epsilon(1);
    end
    %fprintf("HX effectiveness = %f\n", epsilon);
    Qpred = epsilon * min(C_c, C_h) * (T_h_i - T_c_i);
    %fprintf("Qpred = %f\n", Qpred/1e+06);

    
  end

%% Length correction based on Q prediction %%

f = @(L) Q_pred_L(L) - (Q_tot/n_modules); %function that calculates the difference bw predicted and actual heat transfer;
f_scalar = @(L) f(L(1));

L_low = 1e-3;
L_high = 1;  % increases if Q_pred(BD_hi) < Q_tot

% fprintf('f(L_low) = %g\n', f(L_low));
% fprintf('f(L_high) = %g\n', f(L_high));

fprintf("Total heat loss per module = %f\n", Q_tot/n_modules)

while true
    Q_high = Q_pred_L(L_high);
    fprintf("L bracket upper = %.4f m,  Q = %.2f W\n", L_high, Q_high);
    if Q_high >= (Q_tot/n_modules)
        break
    end
    L_high = L_high*2;
    if L_high > 20
        fprintf("Max heat loss possible at L = %.4f m is Q = %.2f W\n", ...
                L_high/2, Q_pred_L(L_high/2));
        error('L_high grew too large — check inputs or UA_unit');
    end
end
  % ylow = f_scalar(L_low);
  % yhigh = f_scalar(L_high);
  % disp([L_low, ylow; L_high, yhigh])

% fprintf('L_low = %g\n', L_low);
% fprintf('L_high) = %g\n', L_high);

%% L_hx solution %%

L_solution = fzero(f_scalar, [L_low, L_high]);% solving for length of hx using a bissection scheme
Q_pred_solution = Q_pred_L(L_solution);
T_cool_out = T_h_i - Q_pred_solution/(m_dot_coolant * cp_EG_50_50(T_mean_h));
fprintf("Effective length of HX is = %f\n",L_solution); %*cosd(hx_theta)
%fprintf("Surface area = %f\n", A_ht_air*L_solution);
fprintf("Predicted coolant outlet temperature = %.4f K\n", T_cool_out);
fprintf("Required coolant outlet temperature  = %.4f K\n", T_h_o);
fprintf("Deviation = %.4f K\n", T_cool_out - T_h_o);
 
 %% Calculation of dP across the HX %%
 
 sr = 2e-06; %m; 2e-06 - OG; surface roughness of aluminium
 f_hx = (1/(-2 * log10(2.7*log10(Re_air)^1.2/Re_air+(sr/d_h_air)/3.71)))^2; %(1/-(2*log10((2.7*(log(Re_air)^1.2))/Re_air)+((sr/d_h_air)/3.71)))^2;
 %fprintf("friction factor = %f\n",f_hx);
 rho_mean_air = 0.5*(rho_air(P3,T_c_i)+rho_air(P3,T_c_o));
 rho_mean_coolant = 0.5*(rho_EG(T_h_i)+rho_EG(T_h_o));
 dp_hx = 0.5*rho_mean_air*v_channel_air*v_channel_air*f_hx*L_solution/d_h_air;

 fprintf("Pressure loss = %f\n",dp_hx);
 
 
 %% Estimation of state variables at end of the HX %%

 A4=A3;
 d4 = d3;

 T4 = T_c_o;
 P4 = P3 - dp_hx;
 fprintf("P4 = %.2f\n",P4);

 if P4 < 0
        error('HX is unfeasible; Check inputs;\n\n');
 end
 %fprintf("rho3 = %f\n",rho_air(P4,T4))
 %fprintf("Mdot3 = %f\n",M_dot_3)
 v4 = M_dot_3/(pi*d3*d3*rho_air(P4,T4)/4);
 a4 = sqrt(gamma_air(cp_air(T4))*R*T4);
 M4 = v4/a4;
 % fprintf("M4 = %f\n",M4);
 % fprintf("v4 = %f\n",v4);
 P4_0 = P4*((1+((gamma_air(cp_air(T4))-1).*M4.^2/2)).^((gamma_air(cp_air(T4)))./(gamma_air(cp_air(T4))-1)));
 T4_0 = T4*(1+((gamma_air(cp_air(T4))-1)*M4.^2)/2);
 M_dot_4 = rho_air(P4,T4)*A4*v4;
 
 F_drag_hx = dp_hx*A_ht_air; %(pi*d4*d4)/4; %assuming air wetted area is the reference area for drag force
 cd = F_drag_hx/(0.5*rho_mean_air*v4*v4);
 %fprintf("P4 = %f\n",P4)
 
 %% Mass calculations %%
 V_plates = (N_pass_tot+2)*t_plate*b_hx*L_solution;
 A_eff_fin_air = (0.5*(b_t_air+(t_fin/sqrt(2)*2))*(h_passage_air)) - (0.5*b_t_air*ht_air);
 A_eff_fin_coolant = (0.5*(b_t_coolant+(t_fin/sqrt(2)*2))*(h_passage_coolant)) - (0.5*b_t_coolant*ht_coolant);
 V_fins_air = ((A_eff_fin_air)*L_solution)*0.5*N_fin_air*N_air_pass;
 V_fins_coolant = ((A_eff_fin_coolant)*b_hx)*0.5*N_fin_air*N_air_pass;
 V_tot = V_plates+V_fins_coolant+V_fins_air;
 M_hx = 2700*V_tot; %total mass of hx. assuming it is made of aluminium hence density of 2700 kg/m3;


% Plates
N_plates = N_pass_tot + 2;    % or +2 for covers 
V_plates = N_plates * t_plate * b_hx * L_solution;


% Air side triangular fins
theta_air = deg2rad(45);     % if your angle is given in degrees

A_xsec_air = b_hx * (t_fin/pitch_fin_air) * (2 * ht_air / sin(theta_air));
A_xsec_cool = b_hx * (t_fin/pitch_fin_coolant) * (2 * ht_coolant / sin(theta_air));

% Volume (extruding along flow length)
V_fins_air = A_xsec_air * L_solution * N_air_pass;
V_fins_cool = A_xsec_cool * L_solution * N_coolant_pass;

V_core  = V_plates + V_fins_air + V_fins_cool ;

% Total
rho_Al = 2700;                % kg/m^3
V_tot  = V_core;
M_hx   = rho_Al * V_tot;


pres_term = (P3*A3) - (P4*A4);
mom_term = M_dot_3*(-v3+v4);
drag_HX = +pres_term-mom_term;
%fprintf("HX drag = %f N \n",drag_HX)

if counter == 1
    P4_for_loop = P4;
end

counter = counter+1;
%fprintf("counter = %f\n",counter);


if counter > 200
    error('exiting HX');
end

f_coolant = (1/(-1.8*log10((0.0015/d_h_coolant)^1.11 + (6.9/Re_coolant))))^2;
dp_coolant_loop = f_coolant*b_hx*(0.5*1082*v_channel_coolant^2)/d_h_coolant;
fprintf("Coolant side pressure loss is %f bar \n",dp_coolant_loop/1e+05);
% if dp_coolant_loop > 1.5e+05
%     error("Coolant side pressure loss is %f bar \n",dp_coolant_loop/1e+05);
% end

A_solid_alt = A_frontal - A_o_air*N_fin_air*N_air_pass;

% fprintf('\nPlate area: %.5f m^2, fin area: %.5f m^2, together: %.5f m^2\n', A_p_air, A_f_air, A_p_air + A_f_air);
% fprintf('Solid area: %.5f m^2\n', A_solid);
% fprintf('Open area: %.5f m^2, together: %.5f m^2\n', A_o_air*N_fin_air*N_air_pass, A_solid+A_o_air*N_fin_air*N_air_pass);
fprintf('Diffuser outlet area: %.5f m^2, HX inlet area: %.5f m^2\n', A3, A_frontal);
fprintf('Cross sectional area of core: %.5f m^2\n', V_core/L_solution)

% - 0.5*rho_air(P4,T4)*v4^2
inlet_dp = (0.5*rho_air(P3,T3)*v3^2)*(A_solid_alt/A_frontal);
fprintf('Estimate inlet bulk pressure loss: %.2f Pa\n\n', inlet_dp)

end


