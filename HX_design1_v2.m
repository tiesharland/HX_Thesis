function [dp_coolant_loop,d_h_air, M_dot_4, b_t_air, b_t_coolant, dp_hx, N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, NTU, R_tot, v_channel_air, v_channel_coolant,d_h_coolant, A_o_coolant, A_o_air, Re_air, Re_coolant, h_air, h_coolant, L_solution, UA_unit, v4, P4_0, M4, T4, T4_0, P4, F_drag_hx, M_hx, A4, drag_HX, Q_pred_solution, T_cool_out, Pr_air, f_air, Nu_air, f_hx, Pr_coolant, f_coolant, Nu_coolant, inlet_dp, outlet_dp] = HX_design1_v2(s_bend_switch,e,r,hx_theta, counter, ~,A3,v3, R,P3,P3_0, h, d3, T_h_o, n_modules, T_h_i, T_c_i, T_c_o, T_mean_h, T_mean_c, Q_tot, C_h, C_c, C_star,Pr_air, Pr_coolant, M_dot_coolant_max, M_dot_3, T3)

%% Inputs %%

A_frontal = pi*d3*d3/(4*cosd(hx_theta)); %m^2; Total frontal area % 0.3969; %
fprintf("HX frontal area = %f\n",A_frontal);
%fprintf("e = %f",e);

m_dot_air = M_dot_3;
%fprintf("m_dot_air in hx = %f\n",m_dot_air)
Q_limit = m_dot_air*cp_air(T3)*(T_h_i-T_c_i);
fprintf("Qlim = %f\n",Q_limit);
if Q_limit < Q_tot/n_modules
   error ("HX mass flow rate not sufficient, increase diffuser inlet area or need a puller fan;");
end
m_dot_coolant = M_dot_coolant_max; %kg/s; Coolant mass flow rate

t_plate = 0.5e-03; %m; thickness of the plate in between consecutive air and coolant passages
t_fin = 0.1e-03; %m; thickness of the fin

b_t_air = 1e-03*e; %1.5e-03; %m; side of the equilateral triangle on air side
b_t_coolant = 1e-03*r; %2e-03; %m; side of the equilateral triangle on coolant side

pitch_fin_air = (0.5*b_t_air)+t_fin; %distance between start of fin and mid point of next fin
pitch_fin_coolant = (0.5*b_t_coolant)+t_fin; %distance between start of fin and mid point of next fin

ht_air = (sqrt(3)/2)*b_t_air; %m; calculation of height of triangular fin
ht_coolant = (sqrt(3)/2)*b_t_coolant;
%fprintf("Height of fin air side;coolant side = %f %f \n",ht_air, ht_coolant);

h_passage_air = ht_air+(t_fin*sqrt(2)); %m; space in between two consecutive plates
h_passage_coolant = ht_coolant+(t_fin*sqrt(2));

%% Calculating N passes %%

b_hx = sqrt(A_frontal); %Assuming the hx frontal area is equal to area of circle at diffuser outlet
h_hx = b_hx;
fprintf("HX frontal area = %f\n",b_hx*h_hx);
N_pass_tot = h_hx/((t_plate+(ht_air+(t_fin/sqrt(2))))+((t_plate+(ht_coolant+(t_fin/sqrt(2)))))); %total number of passages (air+coolant); previous one was without the second term in the denominator - changed: 03/03/26
%fprintf("Passages = %f\n",N_pass_tot);
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


%% Effectiveness calculations %% 


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
epsilon = 0;
f_air = 0; 
Nu_air = 0;
Nu_coolant = 0;
f_coolant = 0;

  function [Qpred] = Q_pred_L(L) 
    
    %----------------- Air side ------------------%  

    A_p_air = (b_t_air*0.5)+((b_t_air*0.5)+(t_fin/sqrt(2))); %m; air side primary heat transfer surface
    A_f_air = sqrt(ht_air^2 + (b_t_air*0.5)^2)+sqrt(((b_t_air*0.5)+(t_fin/sqrt(2)))^2 + (ht_air+(t_fin/sqrt(2)))^2); %m; air side fin surface area
 

    A_ht_air = (((A_f_air*N_fin_air*N_air_pass) + (A_p_air*N_fin_air*N_air_pass))*L); %m^2; total heat transfer area


    A_ht_air_1 = b_t_air*2*N_fin_air*N_air_pass*L; %3*b_t_air*2*N_fin_air*N_air_pass*L;
    
    sr = 2e-06; % m - surface roughness of aluminum
    Re_low  = 1800;
    Re_high = 3000;

        
    % Laminar correlation
    f_lam = ((1.8*log10(Re_air))-1.5)^(-2); %%Konakov, Gnielinski et al as cited in Mortean et.al 2019

    % Turbulent correlation
    f_turb = (1/(-2 * log10( 2.7*log10(Re_air)^1.2/Re_air + (sr/d_h_air)/3.71)))^2; %M.Gerl et.al valid for turbulent

    if Re_air <= Re_low
        f_air = f_lam;

    elseif Re_air >= Re_high
        f_air = f_turb;

    else
        x = (Re_air - Re_low)/(Re_high - Re_low);
        w = 0.5*(1-cos(pi*x));   % cosine smoothing

        f_air = (1-w)*f_lam + w*f_turb;
    end
    K = (T_mean_c/T_mean_h)^0.45;

    if Re_air < 2300
        %Nu_air = 3.66;  % Graetz solution for constant wall temperature
        Nu_air = 4.364 + ((0.086)*(1*Re_air*Pr_air*d_h_air/L)^(1.33))/(1+(0.1*Pr_air*(d_h_air*Re_air/L)^(0.83)));
    elseif Re_air < 3000
        w = (Re_air - 2300)/(3000 - 2300);
        Nu_lam = 4.364 + ((0.086)*(1*Re_air*Pr_air*d_h_air/L)^(1.33))/(1+(0.1*Pr_air*(d_h_air*Re_air/L)^(0.83)));
        Nu_turb = (f_air/8)*(Re_air-1000)*Pr_air*(1+(d_h_air/L)^(2/3))*K/(1+(12.7*sqrt(f_air/8)*(Pr_air^(2/3)-1))); %Gnielinski et.al as cited in Mortean et.al 2019
        Nu_air = (1-w)*Nu_lam + w*Nu_turb;
    else
    % Turbulent Gnielinski
        Nu_air = (f_air/8)*(Re_air-1000)*Pr_air*(1+(d_h_air/L)^(2/3))*K/(1+(12.7*sqrt(f_air/8)*(Pr_air^(2/3)-1))); %Gnielinski et.al as cited in Mortean et.al 2019
    end

    h_air = Nu_air*k_air(T_mean_c)/d_h_air;

    k_fin = 237; %W/mK; thermal conductivity of aluminium
    m_air = sqrt(2*h_air/(k_fin*t_fin)); %constant used in fin efficiency calcs.
    eta_fin_air = tanh(m_air*ht_air)/(m_air*ht_air); %fin efficiency

    R_cond = (t_plate/(k_fin*t_plate*(2*N_pass_tot + 2)*b_hx)); %K/W; conductive heat transfer resistance

    R_air = 1 / (h_air * A_ht_air *eta_fin_air); %K.m/W; convective heat transfer resistance air side

    %----------------- Coolant side ------------------%  

    N_fin_coolant = 2*L/(2*t_fin/sqrt(2)+b_t_coolant); %per passage

    A_p_coolant = (b_t_coolant*0.5)+((b_t_coolant*0.5)+(t_fin/sqrt(2))); %m^2; coolant side primary heat transfer surface
    A_f_coolant = sqrt(ht_coolant^2 + (b_t_coolant*0.5)^2)+sqrt(((b_t_coolant*0.5)+(t_fin/sqrt(2)))^2 + (ht_coolant+(t_fin/sqrt(2)))^2);  %m^2; coolant side fin surface area
    
    A_ht_coolant = (((A_f_coolant*N_fin_coolant*N_coolant_pass) + (A_p_coolant*N_fin_coolant*N_coolant_pass))*b_hx); %m^2; total heat transfer area
    
    v_channel_coolant = (m_dot_pass_coolant/(N_fin_coolant))/(rho_EG(T_mean_h)*(A_o_coolant));
    Re_coolant = rho_EG(T_mean_h)*v_channel_coolant*d_h_coolant/mu_EG(T_mean_h);

    Re_low  = 1800;
    Re_high = 3000;

    % Laminar correlation
    f_lam = ((1.8*log10(Re_coolant))-1.5)^(-2);

    % Turbulent correlation
    f_turb = (1/(-2 * log10( 2.7*log10(Re_coolant)^1.2/Re_coolant + (sr/d_h_coolant)/3.71)))^2;

    if Re_coolant <= Re_low
        f_coolant = f_lam;

    elseif Re_coolant >= Re_high
        f_coolant = f_turb;

    else
        x = (Re_coolant - Re_low)/(Re_high - Re_low);
        w = 0.5*(1-cos(pi*x));   % cosine smoothing

        f_coolant = (1-w)*f_lam + w*f_turb;
    end


    if Re_coolant < 2300
        % ---------------- LAMINAR FLOW ----------------
        %f_coolant = 0.25*((1.8*log10(Re_coolant))-1.5)^(-2); %Konakov, Gnielinski et al as cited in Mortean et.al 2019
        Nu_coolant = 4.364 + ((0.086)*(1*Re_coolant*Pr_coolant*d_h_coolant/b_hx)^(1.33))/(1+(0.1*Pr_coolant*(d_h_coolant*Re_coolant/b_hx)^(0.83))); %Stephan and Preußer cited in Mortean et.al
        %Nu_coolant = 3.66;  % constant wall temperature (equilateral triangular passage)
        % f_coolant = 53.3/Re_coolant;   % laminar, equilateral triangular duct - Darcy
        h_coolant = Nu_coolant * k_EG(T_mean_h) / d_h_coolant;
    
       elseif Re_coolant < 3000
        %sr = 2e-06; %m; 2e-06 - OG; surface roughness of aluminium
        %f_coolant = (1/(-2 * log10(2.7*log10(Re_coolant)^1.2/Re_coolant+(sr/d_h_coolant)/3.71)))^2; %M.Gerl et.al valid for turbulent
        w = (Re_coolant - 2300)/(3000 - 2300);
        K1 = 1; %Pr_mean/Pr_wall; assuming it to be 1 since wall temp is assumed to be bulk temp in calculations;
        Nu_lam_coolant = 4.364 + ((0.086)*(1*Re_coolant*Pr_coolant*d_h_coolant/b_hx)^(1.33))/(1+(0.1*Pr_coolant*(d_h_coolant*Re_coolant/b_hx)^(0.83)));%3.66;
        Nu_turb_coolant = (f_coolant/8)*(Re_coolant-1000)*Pr_coolant*(1+(d_h_coolant/b_hx)^(2/3))*K1/(1+(12.7*sqrt(f_coolant/8)*(Pr_coolant^(2/3)-1))); %Gnielinski et.al as cited in Mortean et.al 2019
        Nu_coolant = (1-w)*Nu_lam_coolant + w*Nu_turb_coolant;
    else
    % Turbulent Gnielinski
        %sr = 2e-06; %m; 2e-06 - OG; surface roughness of aluminium
        K1 = 1; %Pr_mean/Pr_wall; assuming it to be 1 since wall temp is assumed to be bulk temp in calculations;
        %f_coolant = (1/(-2 * log10(2.7*log10(Re_coolant)^1.2/Re_coolant+(sr/d_h_coolant)/3.71)))^2; %M.Gerl et.al valid for turbulent
        Nu_coolant = (f_coolant/8)*(Re_coolant-1000)*Pr_coolant*(1+(d_h_coolant/b_hx)^(2/3))*K1/(1+(12.7*sqrt(f_coolant/8)*(Pr_coolant^(2/3)-1))); %Gnielinski et.al as cited in Mortean et.al 2019

    end
    h_coolant = Nu_coolant * k_EG(T_mean_h) / d_h_coolant;   % set in every regime (was only set when laminar)

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

L_low = 1e-4;
L_high = 1;  % increases if Q_pred(BD_hi) < Q_tot

% fprintf('f(L_low) = %g\n', f(L_low));
% fprintf('f(L_high) = %g\n', f(L_high));

while Q_pred_L(L_high) < (Q_tot/n_modules)
    fprintf("I am here\n");
    L_high = L_high * 2;
    fprintf ("L predicted = %f\n",L_high)
    if L_high > 2
        fprintf("Max heat loss possible is only %f\n",Q_pred_L(L_high));
        error('L_high grew too large — check inputs or UA_unit');
        
    end
end

%   ylow = f_scalar(L_low);
%   yhigh = f_scalar(L_high);
%   disp([L_low, ylow; L_high, yhigh])
% 
% fprintf('L_low = %g\n', L_low);
% fprintf('L_high) = %g\n', L_high);

%% L_hx solution %%

L_solution = fzero(f_scalar, [L_low, L_high]);% solving for length of hx using a bissection scheme
%L_solution = 0.2; %L_solution;
fprintf("Effectvie length of HX is = %f\n",L_solution); %*cosd(hx_theta)
%fprintf("Surface area = %f\n", A_ht_air*L_solution);
Q_pred_solution = Q_pred_L(L_solution);
T_cool_out = T_h_i - Q_pred_solution/(m_dot_coolant * cp_EG_50_50(T_mean_h));
%fprintf("Surface area = %f\n", A_ht_air*L_solution);
fprintf("Predicted coolant outlet temperature = %.4f K\n", T_cool_out);
fprintf("Required coolant outlet temperature  = %.4f K\n", T_h_o);
fprintf("Deviation = %.4f K\n", T_cool_out - T_h_o);
 
 %% Calculation of dP across the HX %%
 
 rho_mean_air = 0.5*(rho_air(P3,T_c_i)+rho_air(P3,T_c_o));
 rho_mean_coolant = 0.5*(rho_EG(T_h_i)+rho_EG(T_h_o));
 dp_fric = 0.5*rho_mean_air*v_channel_air*v_channel_air*f_air*L_solution/d_h_air;

 % Entrance and exit pressure change (Kays & London, 1960), same as HX_design1_disc.
 % inlet_dp is a loss, outlet_dp a pressure rise; dp_hx = friction + inlet - outlet.
 sigma    = A_o_air*0.5*N_fin_air*N_air_pass/A_frontal;   % free-flow / frontal area (0.5 as in v_channel_air)
 G_ch     = rho_air(P3,T_mean_c)*v_channel_air;           % air mass flux in the channels
 v_in     = G_ch/rho_air(P3,T3);
 Kc       = kays_london_triangular(sigma, G_ch*d_h_air/mu_air(T3), 'Kc', false);
 inlet_dp = 0.5*rho_air(P3,T3)*v_in^2*(1 - sigma^2 + Kc);
 T_exit   = T_c_i + epsilon*(T_h_i - T_c_i);              % same as T4 below
 rho_exit = rho_air(P3 - dp_fric - inlet_dp, T_exit);      % exit state estimate
 v_exit   = G_ch/rho_exit;
 Ke       = kays_london_triangular(sigma, G_ch*d_h_air/mu_air(T_exit), 'Ke', false);
 outlet_dp = 0.5*rho_exit*v_exit^2*(1 - sigma^2 - Ke);
 dp_hx = dp_fric + inlet_dp - outlet_dp;
 fprintf('Sigma: %f\n', sigma)
 fprintf('Entrance pressure drop: %.2f Pa\n', inlet_dp)
 fprintf('Exit pressure rise: %.2f Pa\n', outlet_dp)
 fprintf("Pressure loss = %f\n",dp_hx);
 fprintf("Nusselt air = %f; Re air = %f; f_air = %f\n",Nu_air,Re_air,f_air);
 fprintf("Nusselt cool = %f; Re cool = %f; f_cool = %f\n",Nu_coolant,Re_coolant,f_coolant);
 
 %% Estimation of state variables at end of the HX %%

 A4=A3;
 d4 = d3;
 %fprintf("I am using epsilon = %f\n",epsilon);
 T4 = T_c_i+(epsilon*(C_c*(T_h_i-T_c_i)/C_c)); %T_c_o; %re-update with calculated effectiveness instead of using previously used T_c_o
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
 
 if s_bend_switch == "YES"
    K_sbend = 0.6; %Idelchik correlation for S bend chapter 5 pg 442
    dp_S_bend_annular = K_sbend * 0.5 * rho_air(P4,T4) * v4^2;
    P4 = P4 - (2*dp_S_bend_annular);
 end
 
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
theta_air = deg2rad(45);     

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

%f_coolant = (1/(-1.8*log10((0.0015/d_h_coolant)^1.11 + (6.9/Re_coolant))))^2;
sigma_coolant = A_o_coolant*N_fin_coolant*N_coolant_pass/(L_solution*h_hx);
G_coolant = m_dot_coolant/(A_o_coolant*N_fin_coolant*N_coolant_pass);
dp_minor = (0.5*(1-sigma_coolant) + (1-sigma_coolant)^2)*G_coolant*G_coolant/(2*rho_mean_coolant);
fprintf("Minor coolant dP = %f\n",dp_minor);
dp_coolant_loop = f_coolant*b_hx*(0.5*rho_mean_coolant*v_channel_coolant^2)/d_h_coolant;
fprintf("Coolant side pressure loss is %f bar; friction facto = %f \n",dp_coolant_loop/1e+05,f_coolant);

%% Extra outputs, same meaning as in HX_design1 (used by compare_models) %%
f_hx = f_air;   % friction factor used for dp_hx (here the same as f_air)
f_air = f_air/4; % Match outputs -> fanning / skin friction
% inlet_dp and outlet_dp are the Kays & London values, set with dp_hx above
% if dp_coolant_loop > 1.5e+05
%     error("Coolant side pressure loss is %f bar \n",dp_coolant_loop/1e+05);
% end

end