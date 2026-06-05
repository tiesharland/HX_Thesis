function [dp_coolant_loop, d_h_air, M_dot_4, b_t_air, b_t_coolant, dp_hx, ...
          N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, NTU, R_tot, ...
          v_channel_air, v_channel_coolant, d_h_coolant, A_o_coolant, A_o_air, ...
          Re_air, Re_coolant, h_air, h_coolant, L_solution, UA_unit, v4, P4_0, ...
          M4, T4, T4_0, P4, F_drag_hx, M_hx, A4, drag_HX, T_cool_out_arr, dp_cool_seg, ...
          T_air_seg, P_air_seg, v_air_seg, Re_air_seg, Pr_air_seg, v_channel_seg, ...
          K_seg, f_air_seg, Nu_seg_arr, h_air_seg, eta_fin_seg, ...
          UA_seg_arr, NTU_seg_arr, eps_seg_arr, Q_seg_arr, dp_seg_arr, ...
          Q_pred_solution, T_h_o_solution] = ...
          HX_design1_disc(e, r, hx_theta, counter, A3, v3, R, P3, ...
          d3, T_h_o, n_modules, T_h_i, T_c_i, T_mean_h, ...
          Q_tot, M_dot_coolant, M_dot_3, T3, N_segments, tol_T)


% T_c_i = T3;
% T_h_o = T_h_o - 2;

%% ---- Geometry and passage counts ---- %%

sr  = 2e-06;
k_fin_val = 237;

A_frontal = pi*d3*d3/(4*cosd(hx_theta));
fprintf("HX frontal area = %f\n", A_frontal);
fprintf("e = %f\n", e);

m_dot_air = M_dot_3;
Q_limit   = m_dot_air*cp_air(T3)*(T_h_i-T3);
fprintf("Qlim = %f\n", Q_limit);
if Q_limit < Q_tot/n_modules
    error("HX mass flow rate not sufficient, increase diffuser inlet area or need a puller fan;");
end

m_dot_cool_seg = M_dot_coolant / N_segments;

t_plate     = 0.5e-03;
t_fin       = 0.1e-03;
b_t_air     = 1e-03*e;
b_t_coolant = 1e-03*r;

pitch_fin_air     = (0.5*b_t_air)     + t_fin;
pitch_fin_coolant = (0.5*b_t_coolant) + t_fin;

ht_air     = (sqrt(3)/2)*b_t_air;
ht_coolant = (sqrt(3)/2)*b_t_coolant;

h_passage_air     = ht_air     + (t_fin*sqrt(2));
h_passage_coolant = ht_coolant + (t_fin*sqrt(2));

b_hx = sqrt(A_frontal);
h_hx = b_hx;
fprintf("HX height= %f\n", b_hx);

N_pass_tot = h_hx / ((t_plate+(ht_air+(t_fin/sqrt(2)))) + ...
                      (t_plate+(ht_coolant+(t_fin/sqrt(2)))));
fprintf("Passages = %f\n", N_pass_tot);

N_fin_air      = 2*floor(b_hx/((2*t_fin/sqrt(2))+b_t_air));
N_air_pass     = floor(N_pass_tot)+1;
N_coolant_pass = N_air_pass-1;

d_h_air     = 4*(0.5*b_t_air*ht_air)        /(3*b_t_air);
d_h_coolant = 4*(0.5*b_t_coolant*ht_coolant)/(3*b_t_coolant);

A_o_air     = (0.5*(b_t_air*0.5)*ht_air) + ...
              (0.5*((b_t_air*0.5)+(t_fin/sqrt(2)))*(ht_air+(t_fin/sqrt(2))));
A_o_coolant = (0.5*(b_t_coolant*0.5)*ht_coolant) + ...
              (0.5*((b_t_coolant*0.5)+(t_fin/sqrt(2)))*(ht_coolant+(t_fin/sqrt(2))));

m_dot_pass_air     = m_dot_air     / N_air_pass;
m_dot_pass_cool_seg = m_dot_cool_seg / N_coolant_pass;


A_p_coolant   = (b_t_coolant*0.5) + ((b_t_coolant*0.5)+(t_fin/sqrt(2)));
A_f_coolant   = sqrt(ht_coolant^2+(b_t_coolant*0.5)^2) + ...
    sqrt(((b_t_coolant*0.5)+(t_fin/sqrt(2)))^2 + ...
    (ht_coolant+(t_fin/sqrt(2)))^2);

A_p_air      = (b_t_air*0.5) + ((b_t_air*0.5)+(t_fin/sqrt(2)));
A_f_air      = sqrt(ht_air^2+(b_t_air*0.5)^2) + ...
    sqrt(((b_t_air*0.5)+(t_fin/sqrt(2)))^2 + ...
    (ht_air+(t_fin/sqrt(2)))^2);

R_cond_seg = t_plate/(k_fin_val*t_plate*(2*N_pass_tot+2)*b_hx);


%% ---- Inlet bulk reference (scalar outputs) ---- %%

v_channel_air = (m_dot_pass_air/N_fin_air) / (rho_air(P3,T3)*A_o_air);
Re_air        = rho_air(P3,T3)*v_channel_air*d_h_air / mu_air(T3);

%% ---- Shared variables written by nested function ---- %%

T_c_o_disc        = T3;
dp_hx_disc        = 0;
v_channel_coolant = 0;
Re_coolant        = 0;
h_air             = 0;
h_coolant         = 0;
NTU               = 0;
R_tot             = 0;
UA_unit           = 0;

% Node arrays: N_segments+1 points
T_air_seg    = zeros(1, N_segments+1);
P_air_seg    = zeros(1, N_segments+1);
Re_air_seg   = zeros(1, N_segments+1);
Pr_air_seg   = zeros(1, N_segments+1);
f_air_seg    = zeros(1, N_segments+1);
v_channel_seg    = zeros(1, N_segments+1);
v_air_seg    = zeros(1, N_segments+1);
rho_air_seg  = zeros(1, N_segments+1);

% Segment arrays: N_segments values
K_seg        = zeros(1, N_segments);
Q_seg_arr    = zeros(1, N_segments);
dp_seg_arr   = zeros(1, N_segments);
Nu_seg_arr   = zeros(1, N_segments);
h_air_seg    = zeros(1, N_segments);
eta_fin_seg  = zeros(1, N_segments);
UA_seg_arr   = zeros(1, N_segments);
NTU_seg_arr  = zeros(1, N_segments);
eps_seg_arr  = zeros(1, N_segments);
T_cool_out_arr = zeros(1, N_segments);
v_cool_seg   = zeros(1, N_segments);
dp_cool_seg  = zeros(1, N_segments);

%% ---- Nested: counter-flow effectiveness ---- %%

function eps = epsilon_counterflow(NTU, C_star)
    if C_star < 1e-3
        eps = 1 - exp(-NTU);
    else
        eps = (1 - exp(-NTU.*(1-C_star))) ./ ...
              (1 - C_star.*exp(-NTU.*(1-C_star)));
    end
end

%% ---- Nested: discretised Q prediction ---- %%

function [Q_total, T_h_o_predicted] = Q_pred_L_disc(L)
    % Air flows along L, discretised into N_segments equal slices of dx.
    % Coolant flows along b_hx perpendicular to air (cross-flow geometry).
    % All fluid properties evaluated locally at each node temperature.
    % Node arrays track state at N+1 nodes (inlet + segment outlets).
    % Segment arrays track integrated quantities over each segment.

    dx  = L / N_segments;

    Q_total  = 0;
    dp_total = 0;

    tol_cool = 1e-4; % K
    max_cool = 50;

    %% ---- Initialise air inlet node and coolant ---- %%

    T_cool_out = T3;
    rho_in = rho_air(P3, T3);
    v_channel_in = v3*A_frontal / (A_o_air * N_fin_air * N_air_pass)
    v_channel_in = v3*A3 / (A_o_air * N_fin_air * N_air_pass)
    v_channel_in = (m_dot_pass_air/N_fin_air) / (rho_in*A_o_air)
    % v_channel_in = v3;
    A_o_air * N_fin_air * N_air_pass
    A3
    Re_in = rho_air(P3, T3)*v_channel_in*d_h_air / mu_air(T3);
    Pr_in = mu_air(T3)*cp_air(T3)/k_air(T3);
    f_in = (1/(-2*log10(2.7*log10(Re_in)^1.2/Re_in+(sr/d_h_air)/3.71)))^2;
    T_in = T3;
    P_in = P3;
    v_air_in = v3

    T_air_seg(1)   = T3;      P_air_seg(1)   = P3;
    v_channel_seg(1)   = v_channel_in;      Re_air_seg(1)  = Re_in;
    Pr_air_seg(1)  = Pr_in;     f_air_seg(1)   = f_in;     
    v_air_seg(1)   = v3;        rho_air_seg(1) = rho_in;


    %% ---- Geometry ---- %%

    N_fin_cool  = 2*dx/(2*t_fin/sqrt(2)+b_t_coolant);

    A_ht_air_seg = ((A_f_air*N_fin_air*N_air_pass) + ...
        (A_p_air*N_fin_air*N_air_pass))*dx;

    % %% ---- Nested: compute all node properties from (T, P) ---- %%
    % 
    % function node = compute_node(T, P)
    %     node.T   = T;
    %     node.P   = P;
    %     node.rho = rho_air(P, T);
    %     node.v_channel   = (m_dot_pass_air/N_fin_air) / (node.rho*A_o_air);
    %     node.Re  = node.rho*node.v_channel*d_h_air / mu_air(T);
    %     node.Pr  = mu_air(T)*cp_air(T)/k_air(T);
    %     node.f   = (1/(-2*log10(2.7*log10(node.Re)^1.2/node.Re + ...
    %                (sr/d_h_air)/3.71)))^2;
    % end

    %% ---- Segment loop (coolant iteration)---- %%

    for k = 1:N_segments

        for inner = 1:max_cool

            T_mean_h_seg = (T_h_i + T_cool_out) / 2;

            % Coolant properties at T_mean_h_seg

            v_cool      = (m_dot_pass_cool_seg/(2*N_fin_cool)) / (rho_EG(T_mean_h_seg)*A_o_coolant);
            Re_cool     = rho_EG(T_mean_h_seg)*v_cool*d_h_coolant / mu_EG(T_mean_h_seg);
            Pr_cool     = mu_EG(T_mean_h_seg)*cp_EG_50_50(T_mean_h_seg)/k_EG(T_mean_h_seg);
            f_cool      = (1/(-1.8*log10((0.0015/d_h_coolant)^1.11 + (6.9/Re_cool))))^2;
            dp_cool     = f_cool*b_hx*(0.5*1082*v_cool^2)/d_h_coolant;

            if Re_cool < 2300
                Nu_cool    = 4.36;
                h_cool = Nu_cool*k_EG(T_mean_h_seg)/d_h_coolant;
            else
                Nu_cool = 0.023*Re_cool^0.8*Pr_cool^0.4;
                j_cool = Nu_cool/(Re_cool*Pr_cool^(1/3));
                G_cool = rho_EG(T_mean_h_seg)*v_cool;
                h_cool = j_cool*G_cool* ...
                    cp_EG_50_50(T_mean_h_seg)/Pr_cool^(2/3);
            end

            m_cool       = sqrt(2*h_cool/(k_fin_val*t_fin));
            eta_fin_cool = tanh(m_cool*ht_coolant)/(m_cool*ht_coolant);

            A_ht_cool_seg = (A_f_coolant*N_fin_cool*N_coolant_pass) + ...
                (A_p_coolant*N_fin_cool*N_coolant_pass);

            R_cool_seg = 1/(h_cool*A_ht_cool_seg*eta_fin_cool);


            % Segment capacities using T_mean_h_seg

            K       = T_mean_h_seg / T_in;
            C_c     = m_dot_air*cp_air(T_in);
            C_h     = m_dot_cool_seg*cp_EG_50_50(T_mean_h_seg);
            C_min   = min(C_c, C_h);
            C_star  = C_min/max(C_c, C_h);


            % Air-side heat transfer 

            if Re_in < 2300
                Nu_seg = 3.66;
            elseif Re_in < 3000
                w           = (Re_in-2300)/(3000-2300);
                Nu_turb     = (f_in/2)*(Re_in-1000)*Pr_in* ...
                              (1+(d_h_air/(dx*k))^(2/3))*K / ...
                              (1+12.7*sqrt(f_in/2)*(Pr_in^(2/3)-1));
                Nu_seg = (1-w)*3.66 + w*Nu_turb;
            else
                Nu_seg     = (f_in/2)*(Re_in-1000)*Pr_in* ...
                              (1+(d_h_air/(dx*k))^(2/3))*K / ...
                              (1+12.7*sqrt(f_in/2)*(Pr_in^(2/3)-1));
            end

            h_air   = Nu_seg*k_air(T_in)/d_h_air;
            m_air   = sqrt(2*h_air/(k_fin_val*t_fin));
            eta_fin = tanh(m_air*ht_air)/(m_air*ht_air);
            R_air   = 1/(h_air*A_ht_air_seg*eta_fin);

            % Total resistance and heat transfer
            R_tot   = R_air + R_cond_seg + R_cool_seg;
            UA_unit = 1/R_tot;
            NTU     = UA_unit/C_min;
            eps     = epsilon_counterflow(NTU, C_star);

            Q_seg   = eps*C_min*(T_h_i - T_in);
            dp_seg  = f_in*(dx/d_h_air)*0.5*rho_in*v_channel_in^2;
            % dp_seg  = f_in*(dx/d_h_air)*0.5*rho_in*v_air_in^2;

            % Outlet node conditions
            T_out = T_in + Q_seg/(m_dot_air*cp_air(T_in));
            P_out = P_in - dp_seg;
            T_cool_out_old = T_cool_out;
            T_cool_out = T_h_i - Q_seg/(m_dot_cool_seg * cp_EG_50_50(T_mean_h_seg));

            v_air_out = m_dot_air/(pi*d3*d3*rho_air(P_out,T_out)/4);

            rho_out = rho_air(P_out, T_out);
            v_channel_out = (m_dot_pass_air/N_fin_air) / (rho_out*A_o_air);
            Re_out = rho_out*v_channel_out*d_h_air / mu_air(T_out);
            % Re_out = rho_out*v_air_out*d_h_air / mu_air(T_out);
            Pr_out = mu_air(T_out)*cp_air(T_out)/k_air(T_out);
            f_out = (1/(-2*log10(2.7*log10(Re_out)^1.2/Re_out+(sr/d_h_air)/3.71)))^2;


            difference = T_cool_out - T_cool_out_old;
            if abs(difference) >= tol_cool
                if inner == max_cool
                    fprintf("Segment %.0f did not converge, T_cool_out: %.2f K, diff: %.5f  ", k, T_cool_out, difference)
                    break
                end
            else
                % fprintf("Segment %.0f converged in %.0f, T_cool_out: %.2f K   ", k, inner, T_cool_out)
                break
            end

        end

        Q_seg_arr(k)   = Q_seg;
        dp_seg_arr(k)  = dp_seg;
        Nu_seg_arr(k)  = Nu_seg;
        h_air_seg(k)   = h_air;
        eta_fin_seg(k) = eta_fin;
        UA_seg_arr(k)  = UA_unit;
        eps_seg_arr(k) = eps;
        NTU_seg_arr(k) = NTU;
        K_seg(k)       = K;
        T_cool_out_arr(k) = T_cool_out;
        v_cool_seg(k)  = v_cool;
        dp_cool_seg(k)  = dp_cool;

        % Compute and store outlet node at index k+1

        T_air_seg(k+1)   = T_out;      P_air_seg(k+1)   = P_out;
        v_channel_seg(k+1)   = v_channel_out;      Re_air_seg(k+1)  = Re_out;
        Pr_air_seg(k+1)  = Pr_out;     f_air_seg(k+1)   = f_out;
        v_air_seg(k+1)   = v_air_out;  rho_air_seg(k+1) = rho_out;

        % Update shared scalar outputs
        v_channel_coolant = v_cool;
        Re_coolant        = Re_cool;
        h_coolant         = h_cool;

        T_in   = T_out;
        P_in   = P_out;
        rho_in = rho_out;
        v_channel_in = v_channel_out;
        Re_in = Re_out;
        Pr_in = Pr_out;
        f_in  = f_out;

        Q_total  = Q_total  + Q_seg;
        dp_total = dp_total + dp_seg;

    end


% %% ---- Segment loop (approximate next T_cool_out)---- %%
% 
%     for k = 1:N_segments
% 
%         T_mean_h_seg = (T_h_i + T_cool_out) / 2;
% 
%         % Coolant properties at T_mean_h_seg
% 
%         v_cool      = (m_dot_pass_cool_seg/(2*N_fin_cool)) / (rho_EG(T_mean_h_seg)*A_o_coolant);
%         Re_cool     = rho_EG(T_mean_h_seg)*v_cool*d_h_coolant / mu_EG(T_mean_h_seg);
%         Pr_cool     = mu_EG(T_mean_h_seg)*cp_EG_50_50(T_mean_h_seg)/k_EG(T_mean_h_seg);
%         f_cool      = (1/(-1.8*log10((0.0015/d_h_coolant)^1.11 + (6.9/Re_cool))))^2;
%         dp_cool     = f_cool*b_hx*(0.5*1082*v_cool^2)/d_h_coolant;
% 
%         if Re_cool < 2300
%             Nu_cool    = 4.36;
%             h_cool = Nu_cool*k_EG(T_mean_h_seg)/d_h_coolant;
%         else
%             Nu_cool = 0.023*Re_cool^0.8*Pr_cool^0.4;
%             j_cool = Nu_cool/(Re_cool*Pr_cool^(1/3));
%             G_cool = rho_EG(T_mean_h_seg)*v_cool;
%             h_cool = j_cool*G_cool* ...
%                 cp_EG_50_50(T_mean_h_seg)/Pr_cool^(2/3);
%         end
% 
%         m_cool       = sqrt(2*h_cool/(k_fin_val*t_fin));
%         eta_fin_cool = tanh(m_cool*ht_coolant)/(m_cool*ht_coolant);
% 
%         A_ht_cool_seg = (A_f_coolant*N_fin_cool*N_coolant_pass) + ...
%             (A_p_coolant*N_fin_cool*N_coolant_pass);
% 
%         R_cool_seg = 1/(h_cool*A_ht_cool_seg*eta_fin_cool);
% 
% 
%         % Segment capacities using T_mean_h_seg
% 
%         K       = T_mean_h_seg / T_in;
%         C_c     = m_dot_air*cp_air(T_in);
%         C_h     = m_dot_cool_seg*cp_EG_50_50(T_mean_h_seg);
%         C_min   = min(C_c, C_h);
%         C_star  = C_min/max(C_c, C_h);
% 
% 
%         % Air-side heat transfer 
% 
%         if Re_in < 2300
%             Nu_seg = 3.66;
%         elseif Re_in < 3000
%             w           = (Re_in-2300)/(3000-2300);
%             Nu_turb     = (f_in/2)*(Re_in-1000)*Pr_in* ...
%                           (1+(d_h_air/(dx*k))^(2/3))*K / ...
%                           (1+12.7*sqrt(f_in/2)*(Pr_in^(2/3)-1));
%             Nu_seg = (1-w)*3.66 + w*Nu_turb;
%         else
%             Nu_seg     = (f_in/2)*(Re_in-1000)*Pr_in* ...
%                           (1+(d_h_air/(dx*k))^(2/3))*K / ...
%                           (1+12.7*sqrt(f_in/2)*(Pr_in^(2/3)-1));
%         end
% 
%         h_air   = Nu_seg*k_air(T_in)/d_h_air;
%         m_air   = sqrt(2*h_air/(k_fin_val*t_fin));
%         eta_fin = tanh(m_air*ht_air)/(m_air*ht_air);
%         R_air   = 1/(h_air*A_ht_air_seg*eta_fin);
% 
%         % Total resistance and heat transfer
%         R_tot   = R_air + R_cond_seg + R_cool_seg;
%         UA_unit = 1/R_tot;
%         NTU     = UA_unit/C_min;
%         eps     = epsilon_counterflow(NTU, C_star);
% 
%         Q_seg   = eps*C_min*(T_h_i - T_in);
%         % dp_seg  = f_in*(dx/d_h_air)*0.5*rho_in*v_channel_in^2;
%         dp_seg  = f_in*(dx/d_h_air)*0.5*rho_in*v_air_in^2;
% 
%         % Outlet node conditions
%         T_out = T_in + Q_seg/(m_dot_air*cp_air(T_in));
%         P_out = P_in - dp_seg;
%         T_cool_out_old = T_cool_out;
%         T_cool_out = T_h_i - Q_seg/(m_dot_cool_seg * cp_EG_50_50(T_mean_h_seg));
% 
%         v_air_out = m_dot_air/(pi*d3*d3*rho_air(P_out,T_out)/4);
% 
%         rho_out = rho_air(P_out, T_out);
%         v_channel_out = (m_dot_pass_air/N_fin_air) / (rho_out*A_o_air);
%         % Re_out = rho_out*v_channel_out*d_h_air / mu_air(T_out);
%         Re_out = rho_out*v_air_out*d_h_air / mu_air(T_out);
%         Pr_out = mu_air(T_out)*cp_air(T_out)/k_air(T_out);
%         f_out = (1/(-2*log10(2.7*log10(Re_out)^1.2/Re_out+(sr/d_h_air)/3.71)))^2;
% 
% 
%         % difference = T_cool_out - T_cool_out_old;
%         % if abs(difference) >= tol_cool
%         %     if inner == max_cool
%         %         fprintf("Segment %.0f did not converge, T_cool_out: %.2f K, diff: %.5f  ", k, T_cool_out, difference)
%         %         break
%         %     end
%         % else
%         %     % fprintf("Segment %.0f converged in %.0f, T_cool_out: %.2f K   ", k, inner, T_cool_out)
%         %     break
%         % end
% 
%         Q_seg_arr(k)   = Q_seg;
%         dp_seg_arr(k)  = dp_seg;
%         Nu_seg_arr(k)  = Nu_seg;
%         h_air_seg(k)   = h_air;
%         eta_fin_seg(k) = eta_fin;
%         UA_seg_arr(k)  = UA_unit;
%         eps_seg_arr(k) = eps;
%         NTU_seg_arr(k) = NTU;
%         K_seg(k)       = K;
%         T_cool_out_arr(k) = T_cool_out;
%         v_cool_seg(k)  = v_cool;
%         dp_cool_seg(k)  = dp_cool;
% 
%         % Compute and store outlet node at index k+1
% 
%         T_air_seg(k+1)   = T_out;      P_air_seg(k+1)   = P_out;
%         v_channel_seg(k+1)   = v_channel_out;      Re_air_seg(k+1)  = Re_out;
%         Pr_air_seg(k+1)  = Pr_out;     f_air_seg(k+1)   = f_out;
%         v_air_seg(k+1)   = v_air_out;  rho_air_seg(k+1) = rho_out;
% 
%         % Update shared scalar outputs
%         v_channel_coolant = v_cool;
%         Re_coolant        = Re_cool;
%         h_coolant         = h_cool;
% 
%         T_in   = T_out;
%         P_in   = P_out;
%         rho_in = rho_out;
%         v_channel_in = v_channel_out;
%         Re_in = Re_out;
%         Pr_in = Pr_out;
%         f_in  = f_out;
% 
%         Q_total  = Q_total  + Q_seg;
%         dp_total = dp_total + dp_seg;
% 
%     end

    T_c_o_disc = T_air_seg(end);
    dp_hx_disc = dp_total;

    fprintf("\nQ_total = %.2f W,  T_air_out = %.2f K\n", Q_total, T_air_seg(end));

    % T_h_o_mix = T_h_i - Q_total/(M_dot_coolant*cp_EG_50_50(T_mean_h));
    % fprintf("Mixed coolant outlet = %.2f K,  system T_h_o = %.2f K,  deviation = %.4f K\n", ...
    %         T_h_o_mix, T_h_o, T_h_o_mix - T_h_o);

    % Mixed coolant outlet temperature
    % T_h_o_predicted = sum(T_cool_out_arr .* arrayfun(@cp_EG_50_50, T_cool_out_arr)) / ...
    %     sum(arrayfun(@cp_EG_50_50, T_cool_out_arr));

    T_h_o_predicted = mean(T_cool_out_arr);

    fprintf("Predicted mixed coolant outlet = %.4f K,  T_h_o requirement = %.4f K,  deviation = %.4f K\n", ...
        T_h_o_predicted, T_h_o, T_h_o_predicted - T_h_o);

end

function Q = Q_only(L)
    [Q, ~] = Q_pred_L_disc(L);
end

function T = T_only(L)
    [~, T] = Q_pred_L_disc(L);
end

%% ---- Length iteration ---- %%

fprintf("Total heat loss per module = %f W\n", Q_tot/n_modules);

% tol_T = 7;
L_low  = 1e-3;
L_high = 1;

% Q_satisfied = false;
% T_satisfied = false;

while true

    [Q_high, T_high] = Q_pred_L_disc(L_high);
    fprintf("L bracket upper = %.4f m,  Q = %.2f W,  T_h_o_pred = %.4f K\n", ...
        L_high, Q_high, T_high);

    Q_satisfied = Q_high >= Q_tot/n_modules;
    T_satisfied = T_high <= T_h_o;
    T_close = abs(T_high - T_h_o) <= tol_T;

    if T_satisfied && Q_satisfied
        % Case 1: both satisfied — prefer T_h_o residual
        fprintf("Case 1: both constraints satisfied at L_high = %.4f m\n", L_high);
        break

    elseif Q_satisfied && ~T_satisfied
        % Case 2: Q satisfied but T_h_o not reached
        if T_close
            fprintf('Case 2a: Q satisfied and T within tol_T at L_high = %.4f m\n', L_high);
            break
        else
            fprintf('Case 2b: Q satisfied and T not within tol_T — increasing L_high\n');
        end
        % fprintf("Case 2: Q satisfied but T_h_o not reached — solving on Q residual\n");

    elseif T_satisfied && ~Q_satisfied
        % Case 3: T satisfied but Q not — increase L until Q is also satisfied
        fprintf("Case 3: T_h_o satisfied but Q insufficient — increasing L_high\n");

    else
        % Case 4: neither satisfied — increase L_high
        fprintf("Case 4: neither constraint satisfied — increasing L_high\n");

    end

    L_high = L_high*2;
    if L_high > 20
        [Q_best, T_best] = Q_pred_L_disc(L_high/2);
        fprintf("Max achievable: Q = %.2f W, T_h_o_pred = %.4f K at L = %.4f m\n", ...
            Q_best, T_best, L_high/2);
        error('L_high grew too large — check inputs or UA_unit');
    end

end

if T_satisfied && Q_satisfied
    % Case 1: T_h_o constraint met — solve on T_h_o residual
    % options    = optimset('TolFun', tol_L_iter);
    T_residual = @(L) T_only(L) - T_h_o;
    L_solution = fzero(T_residual, [L_low, L_high]);
    fprintf("L solved on T_h_o residual: L = %.4f m\n", L_solution);
elseif T_close && Q_satisfied
    % Case 2a: T within tolerance of T_h_o - solve on Q residual
    % T_residual = @(L) T_only(L) - (T_h_o + tol_T);
    % L_solution = fzero(T_residual, [L_low, L_high]);
    % fprintf("L solved on T_h_o + %.1f K residual: L = %.4f m\n", tol_T, L_solution);
    Q_residual = @(L) Q_only(L) - Q_tot/n_modules;
    L_solution = fzero(Q_residual, [L_low, L_high]);
    fprintf("L solved on Q residual: L = %.4f m\n", L_solution);
elseif Q_satisfied
    % Case 2b: T_h_o not reachable — solve on Q residual
    Q_residual = @(L) Q_only(L) - Q_tot/n_modules;
    L_solution = fzero(Q_residual, [L_low, L_high]);
    fprintf("L solved on Q residual: L = %.4f m\n", L_solution);
end

fprintf("Effective length of HX = %.4f m\n", L_solution);
[Q_pred_solution, T_h_o_solution] = Q_pred_L_disc(L_solution);

% while true
%     Q_high = Q_pred_L_disc(L_high);
%     fprintf("L bracket upper = %.4f m,  Q = %.2f W\n", L_high, Q_high);
%     if Q_high >= (Q_tot/n_modules)
%         break
%     end
%     L_high = L_high*2;
%     if L_high > 20
%         fprintf("Max heat loss possible at L = %.4f m is Q = %.2f W\n", ...
%                 L_high/2, Q_pred_L_disc(L_high/2));
%         error('L_high grew too large — check inputs or UA_unit');
%     end
% end
% 
% L_solution = fzero(@(L) Q_pred_L_disc(L) - (Q_tot/n_modules), [L_low, L_high]);

% while true
%     T_out_h_high = Q_pred_L_disc(L_high);
%     fprintf("L bracket upper = %.4f m,  T = %.2f K", L_high, T_out_h_high);
%     if T_out_h_high <= T_h_o
%         break
%     end
%     L_high = L_high*2;
%     if L_high > 4
%         fprintf("Max outlet coolant temp possible at L = %.4f m is T = %.2f K", ...
%             L_high/2, Q_pred_L_disc(L_high/2));
%         error('L_high grew too large — check inputs or UA_unit');
%     end
% end
% 
% options = optimset('TolFun', tol_L_iter);
% L_solution = fzero(@(L) Q_pred_L_disc(L) - (T_h_o), [L_low, L_high], options);


% 
% fprintf("Effective length of HX = %.4f m\n", L_solution);
% 
% Q_pred_L_disc(L_solution);

dp_hx = dp_hx_disc;
fprintf("Total pressure loss = %.2f Pa\n", dp_hx);

%% ---- Exit state variables ---- %%

A4  = A3;
T4  = T_c_o_disc;
P4  = P3 - dp_hx;
fprintf("P4 = %.2f Pa\n", P4);

if P4 < 0
    error('HX is unfeasible; Check inputs');
end

v4      = m_dot_air/(pi*d3*d3*rho_air(P4,T4)/4)
v_channel_seg(end)
a4      = sqrt(gamma_air(cp_air(T4))*R*T4);
M4      = v4/a4;
P4_0    = P4*((1+((gamma_air(cp_air(T4))-1).*M4.^2/2)).^ ...
              ((gamma_air(cp_air(T4)))./(gamma_air(cp_air(T4))-1)));
T4_0    = T4*(1+((gamma_air(cp_air(T4))-1)*M4.^2)/2);
M_dot_4 = rho_air(P4,T4)*A4*v4;

%% ---- HX drag ---- %%

A_p_air_out    = (b_t_air*0.5) + ((b_t_air*0.5)+(t_fin/sqrt(2)));
A_f_air_out    = sqrt(ht_air^2+(b_t_air*0.5)^2) + ...
                 sqrt(((b_t_air*0.5)+(t_fin/sqrt(2)))^2 + ...
                      (ht_air+(t_fin/sqrt(2)))^2);
A_ht_air_total = ((A_f_air_out*N_fin_air*N_air_pass) + ...
                  (A_p_air_out*N_fin_air*N_air_pass))*L_solution;

F_drag_hx = dp_hx*A_ht_air_total;

pres_term = (P3*A3) - (P4*A4);
mom_term  = m_dot_air*(-v3+v4);
drag_HX   = pres_term - mom_term;

%% ---- Mass of HX ---- %%

N_plates    = N_pass_tot+2;
V_plates    = N_plates*t_plate*b_hx*L_solution;
theta_hx    = deg2rad(45);
A_xsec_air  = b_hx*(t_fin/pitch_fin_air)*(2*ht_air/sin(theta_hx));
A_xsec_cool = b_hx*(t_fin/pitch_fin_coolant)*(2*ht_coolant/sin(theta_hx));
V_fins_air  = A_xsec_air*L_solution*N_air_pass;
V_fins_cool = A_xsec_cool*L_solution*N_coolant_pass;
V_core      = V_plates+V_fins_air+V_fins_cool;
M_hx        = 2700*V_core;

%% ---- Coolant-side pressure drop ---- %%

% counter = counter+1;
% if counter > 200
%     error('exiting HX');
% end

N_fin_coolant   = 2*L_solution/(2*t_fin/sqrt(2)+b_t_coolant);
f_coolant       = (1/(-1.8*log10((0.0015/d_h_coolant)^1.11 + (6.9/Re_coolant))))^2;
% dp_coolant_loop = f_coolant*b_hx*(0.5*1082*v_channel_coolant^2)/d_h_coolant;
dp_coolant_loop = mean(dp_cool_seg);
fprintf("Coolant side pressure loss = %.4f bar\n", dp_coolant_loop/1e5);

end