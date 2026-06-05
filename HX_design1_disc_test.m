function [dp_coolant_loop, d_h_air, M_dot_4, b_t_air, b_t_coolant, dp_hx, ...
          N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, NTU, R_tot, ...
          v_channel_air, v_channel_coolant, d_h_coolant, A_o_coolant, A_o_air, ...
          Re_air, Re_coolant, h_air, h_coolant, L_solution, UA_unit, v4, P4_0, ...
          M4, T4, T4_0, P4, F_drag_hx, M_hx, A4, drag_HX, ...
          T_air_seg, P_air_seg, v_air_seg, Re_air_seg, Pr_air_seg, ...
          K_seg, f_air_seg, Nu_seg_arr, h_air_seg, eta_fin_seg, ...
          UA_seg_arr, NTU_seg_arr, eps_seg_arr, Q_seg_arr, dp_seg_arr, ...
          T_cool_seg] = ...
          HX_design1_disc_test(e, r, hx_theta, counter, A3, v3, R, P3, ...
          d3, T_h_o, n_modules, T_h_i, T_c_i, T_mean_h, ...
          Q_tot, M_dot_coolant_max, M_dot_3, T3, N_segments)

%% ---- Geometry and passage counts ---- %%

A_frontal = pi*d3*d3/(4*cosd(hx_theta));
fprintf("HX frontal area = %f\n", A_frontal);
fprintf("e = %f\n", e);

m_dot_air = M_dot_3;
Q_limit   = m_dot_air*cp_air(T3)*(T_h_i-T_c_i);
fprintf("Qlim = %f\n", Q_limit);
if Q_limit < Q_tot/n_modules
    error("HX mass flow rate not sufficient, increase diffuser inlet area or need a puller fan;");
end

m_dot_coolant = M_dot_coolant_max;

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
fprintf("HX side length = %f\n", b_hx);

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
m_dot_pass_coolant = m_dot_coolant / N_coolant_pass;

%% ---- Inlet bulk reference (scalar outputs) ---- %%

v_channel_air = (m_dot_pass_air/N_fin_air) / (rho_air(P3,T_c_i)*A_o_air);
Re_air        = rho_air(P3,T_c_i)*v_channel_air*d_h_air / mu_air(T_c_i);

%% ---- Shared variables written by nested function ---- %%

T_c_o_disc        = T_c_i;
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
v_air_seg    = zeros(1, N_segments+1);
Re_air_seg   = zeros(1, N_segments+1);
Pr_air_seg   = zeros(1, N_segments+1);
K_seg        = zeros(1, N_segments+1);
f_air_seg    = zeros(1, N_segments+1);
T_cool_seg   = zeros(1, N_segments+1);

% Segment arrays: N_segments values
Q_seg_arr    = zeros(1, N_segments);
dp_seg_arr   = zeros(1, N_segments);
Nu_seg_arr   = zeros(1, N_segments);
h_air_seg    = zeros(1, N_segments);
eta_fin_seg  = zeros(1, N_segments);
UA_seg_arr   = zeros(1, N_segments);
NTU_seg_arr  = zeros(1, N_segments);
eps_seg_arr  = zeros(1, N_segments);

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

function [Q_total] = Q_pred_L_disc(L)
    % Counter-flow geometry:
    % - Air flows from k=1 (inlet) to k=N (outlet) along L
    % - Coolant flows from k=N+1 (inlet, T_h_i) to k=1 (outlet, T_h_o)
    %   in the opposite direction to air
    %
    % T_cool indexing:
    %   T_cool(k+1) = coolant temperature at inlet to segment k
    %                 (from coolant flow perspective, entering from high-k side)
    %   T_cool(k)   = coolant temperature at outlet of segment k
    %   T_cool(N+1) = T_h_i  (known coolant inlet, at air outlet side)
    %   T_cool(1)   = T_h_o  (known coolant outlet, at air inlet side)
    %
    % Inner iteration converges T_cool distribution by alternating:
    %   1. Air forward march  (k=1 to N) using current T_cool
    %   2. Coolant reverse march (k=N to 1) using updated T_air
    % until T_cool converges.

    dx        = L / N_segments;
    sr        = 2e-06;
    k_fin_val = 237;

    Q_total  = 0;
    dp_total = 0;

    %% ---- Fixed geometry quantities ---- %%

    % Air-side heat transfer area per segment
    A_p_air      = (b_t_air*0.5) + ((b_t_air*0.5)+(t_fin/sqrt(2)));
    A_f_air      = sqrt(ht_air^2+(b_t_air*0.5)^2) + ...
                   sqrt(((b_t_air*0.5)+(t_fin/sqrt(2)))^2 + ...
                        (ht_air+(t_fin/sqrt(2)))^2);
    A_ht_air_seg = ((A_f_air*N_fin_air*N_air_pass) + ...
                    (A_p_air*N_fin_air*N_air_pass))*dx;

    % Coolant-side geometry per segment
    % N_fin_cool: coolant fins in one segment of length dx
    N_fin_cool    = 2*dx/(2*t_fin/sqrt(2)+b_t_coolant);
    A_p_coolant   = (b_t_coolant*0.5) + ((b_t_coolant*0.5)+(t_fin/sqrt(2)));
    A_f_coolant   = sqrt(ht_coolant^2+(b_t_coolant*0.5)^2) + ...
                    sqrt(((b_t_coolant*0.5)+(t_fin/sqrt(2)))^2 + ...
                         (ht_coolant+(t_fin/sqrt(2)))^2);
    A_ht_cool_seg = (A_f_coolant*N_fin_cool*N_coolant_pass) + ...
                    (A_p_coolant*N_fin_cool*N_coolant_pass);

    % Conduction resistance — constant
    R_cond_seg = t_plate/(k_fin_val*t_plate*(2*N_pass_tot+2)*b_hx);

    %% ---- Nested: compute air node properties from (T_air, P_air, T_cool) ---- %%
    % T_cool is the local coolant temperature at this node, used for
    % K correction and capacity rates

    function node = compute_node(T_air, P_air, T_cool_local)
        node.T   = T_air;
        node.P   = P_air;
        node.rho = rho_air(P_air, T_air);
        node.v   = (m_dot_pass_air/N_fin_air) / (node.rho*A_o_air);
        node.Re  = node.rho*node.v*d_h_air / mu_air(T_air);
        node.Pr  = mu_air(T_air)*cp_air(T_air)/k_air(T_air);

        % K uses local coolant temperature, not global T_mean_h
        node.K   = T_cool_local / T_air;

        node.f   = (1/(-2*log10(2.7*log10(node.Re)^1.2/node.Re + ...
                   (sr/d_h_air)/3.71)))^2;

        % Capacity rates use local temperatures
        C_c          = M_dot_3*cp_air(T_air);
        C_h          = M_dot_coolant_max*cp_EG_50_50(T_cool_local);
        node.C_min   = min(C_c, C_h);
        node.C_star  = node.C_min/max(C_c, C_h);
    end

    %% ---- Nested: compute coolant node properties from T_cool ---- %%

    function cool = compute_cool_node(T_cool_local)
        cool.T     = T_cool_local;
        cool.v     = (m_dot_pass_coolant/(2*N_fin_cool)) / ...
                     (rho_EG(T_cool_local)*A_o_coolant);
        cool.Re    = rho_EG(T_cool_local)*cool.v*d_h_coolant / ...
                     mu_EG(T_cool_local);
        cool.Pr    = mu_EG(T_cool_local)*cp_EG_50_50(T_cool_local) / ...
                     k_EG(T_cool_local);

        if cool.Re < 2300
            Nu_c    = 4.36;
            cool.h  = Nu_c*k_EG(T_cool_local)/d_h_coolant;
        else
            Nu_c    = 0.023*cool.Re^0.8*cool.Pr^0.4;
            j_c     = Nu_c/(cool.Re*cool.Pr^(1/3));
            G_c     = rho_EG(T_cool_local)*cool.v;
            cool.h  = j_c*G_c*cp_EG_50_50(T_cool_local)/cool.Pr^(2/3);
        end

        m_c         = sqrt(2*cool.h/(k_fin_val*t_fin));
        cool.eta    = tanh(m_c*ht_coolant)/(m_c*ht_coolant);
        cool.R      = 1/(cool.h*A_ht_cool_seg*cool.eta);
    end

    %% ---- Initial guess for coolant temperature distribution ---- %%
    % T_cool(k+1) = coolant temperature at inlet to air segment k
    % Counter-flow: T_cool increases from k=1 to k=N+1
    % T_cool(1)   = T_h_o (coolant outlet, at air inlet side)
    % T_cool(N+1) = T_h_i (coolant inlet, at air outlet side)

    T_cool = linspace(T_h_o, T_h_i, N_segments+1);

    %% ---- Inner iteration: converge T_cool distribution ---- %%

    tol_inner = 1e-3;   % K
    max_iter  = 100;

    for iter = 1:max_iter

        T_cool_prev = T_cool;

        %--- Air forward march (k=1 to N) ---%
        % T_cool(k+1) is the coolant inlet temperature for segment k
        T_air_iter    = zeros(1, N_segments+1);
        P_air_iter    = zeros(1, N_segments+1);
        T_air_iter(1) = T_c_i;
        P_air_iter(1) = P3;
        Q_iter        = zeros(1, N_segments);

        for k = 1:N_segments
            T_cool_k = T_cool(k+1);  % coolant entering segment k
            node_k   = compute_node(T_air_iter(k), P_air_iter(k), T_cool_k);
            cool_k   = compute_cool_node(T_cool_k);

            R_cool_seg_k = cool_k.R;

            % Nu for this segment
            if node_k.Re < 2300
                Nu_k = 3.66;
            elseif node_k.Re < 3000
                w        = (node_k.Re-2300)/(3000-2300);
                Nu_turb  = (node_k.f/2)*(node_k.Re-1000)*node_k.Pr * ...
                            (1+(d_h_air/(dx*k))^(2/3))*node_k.K / ...
                            (1+12.7*sqrt(node_k.f/2)*(node_k.Pr^(2/3)-1));
                Nu_k     = (1-w)*3.66 + w*Nu_turb;
            else
                Nu_k     = (node_k.f/2)*(node_k.Re-1000)*node_k.Pr * ...
                            (1+(d_h_air/(dx*k))^(2/3))*node_k.K / ...
                            (1+12.7*sqrt(node_k.f/2)*(node_k.Pr^(2/3)-1));
            end

            h_air_k   = Nu_k*k_air(node_k.T)/d_h_air;
            m_air_k   = sqrt(2*h_air_k/(k_fin_val*t_fin));
            eta_fin_k = tanh(m_air_k*ht_air)/(m_air_k*ht_air);
            R_air_k   = 1/(h_air_k*A_ht_air_seg*eta_fin_k);
            R_tot_k   = R_air_k + R_cond_seg + R_cool_seg_k;
            UA_k      = 1/R_tot_k;
            NTU_k     = UA_k/node_k.C_min;
            eps_k     = epsilon_counterflow(NTU_k, node_k.C_star);

            % Driving difference: local coolant inlet minus local air inlet
            Q_iter(k) = eps_k*node_k.C_min*(T_cool_k - node_k.T);

            % Air outlet
            T_air_iter(k+1) = node_k.T + Q_iter(k)/(M_dot_3*cp_air(node_k.T));
            dp_k            = node_k.f*(dx/d_h_air)*0.5*node_k.rho*node_k.v^2;
            P_air_iter(k+1) = P_air_iter(k) - dp_k;
        end

        %--- Coolant reverse march (k=N to 1) ---%
        % Coolant enters at T_cool(N+1) = T_h_i (fixed boundary)
        % and cools as it moves toward k=1
        T_cool(N_segments+1) = T_h_i;  % enforce known boundary

        for k = N_segments:-1:1
            % Coolant flows counter to air: from segment N down to segment 1
            % Coolant inlet to segment k is T_cool(k+1)
            % Energy balance on coolant side: coolant loses Q_iter(k)
            T_cool(k) = T_cool(k+1) - Q_iter(k) / ...
                        (M_dot_coolant_max*cp_EG_50_50(T_cool(k+1)));
        end

        % Enforce known coolant outlet boundary — compare against T_h_o
        % as a convergence diagnostic (not enforced as hard constraint
        % since Q_iter drives the distribution)
        fprintf("Inner iter %d: T_cool(1) = %.4f K, T_h_o = %.4f K, residual = %.6f K\n", ...
                iter, T_cool(1), T_h_o, T_cool(1)-T_h_o);

        % Check convergence of T_cool distribution
        if max(abs(T_cool - T_cool_prev)) < tol_inner
            fprintf("Inner coolant iteration converged in %d iterations\n", iter);
            break
        end

        if iter == max_iter
            warning("Inner coolant iteration did not converge at L=%.4f m", L);
        end

    end

    %% ---- Final air march with converged T_cool ---- %%
    % Populate all segment and node arrays

    node = compute_node(T_c_i, P3, T_cool(2));

    T_air_seg(1)  = node.T;   P_air_seg(1)  = node.P;
    v_air_seg(1)  = node.v;   Re_air_seg(1) = node.Re;
    Pr_air_seg(1) = node.Pr;  K_seg(1)      = node.K;
    f_air_seg(1)  = node.f;
    T_cool_seg(1) = T_cool(1);

    for k = 1:N_segments

        T_cool_k = T_cool(k+1);  % coolant inlet to segment k
        cool_k   = compute_cool_node(T_cool_k);

        R_cool_seg_k = cool_k.R;

        if node.Re < 2300
            Nu_seg = 3.66;
        elseif node.Re < 3000
            w       = (node.Re-2300)/(3000-2300);
            Nu_turb = (node.f/2)*(node.Re-1000)*node.Pr * ...
                      (1+(d_h_air/(dx*k))^(2/3))*node.K / ...
                      (1+12.7*sqrt(node.f/2)*(node.Pr^(2/3)-1));
            Nu_seg  = (1-w)*3.66 + w*Nu_turb;
        else
            Nu_seg  = (node.f/2)*(node.Re-1000)*node.Pr * ...
                      (1+(d_h_air/(dx*k))^(2/3))*node.K / ...
                      (1+12.7*sqrt(node.f/2)*(node.Pr^(2/3)-1));
        end

        h_air_k   = Nu_seg*k_air(node.T)/d_h_air;
        m_air_k   = sqrt(2*h_air_k/(k_fin_val*t_fin));
        eta_fin_k = tanh(m_air_k*ht_air)/(m_air_k*ht_air);
        R_air_k   = 1/(h_air_k*A_ht_air_seg*eta_fin_k);
        R_tot_k   = R_air_k + R_cond_seg + R_cool_seg_k;
        UA_k      = 1/R_tot_k;
        NTU_k     = UA_k/node.C_min;
        eps_k     = epsilon_counterflow(NTU_k, node.C_star);

        Q_seg  = eps_k*node.C_min*(T_cool_k - node.T);
        dp_seg = node.f*(dx/d_h_air)*0.5*node.rho*node.v^2;

        % Store segment arrays
        Q_seg_arr(k)   = Q_seg;
        dp_seg_arr(k)  = dp_seg;
        Nu_seg_arr(k)  = Nu_seg;
        h_air_seg(k)   = h_air_k;
        eta_fin_seg(k) = eta_fin_k;
        UA_seg_arr(k)  = UA_k;
        NTU_seg_arr(k) = NTU_k;
        eps_seg_arr(k) = eps_k;

        Q_total  = Q_total  + Q_seg;
        dp_total = dp_total + dp_seg;

        % Outlet node
        T_out = node.T + Q_seg/(M_dot_3*cp_air(node.T));
        P_out = node.P - dp_seg;

        % Use T_cool(k+1) for the outlet node K and C_star
        % (next segment's coolant inlet is T_cool(k+2), but for the
        %  outlet node itself we use T_cool(k+1) as representative)
        if k < N_segments
            node = compute_node(T_out, P_out, T_cool(k+2));
        else
            node = compute_node(T_out, P_out, T_cool(N_segments+1));
        end

        T_air_seg(k+1)  = node.T;   P_air_seg(k+1)  = node.P;
        v_air_seg(k+1)  = node.v;   Re_air_seg(k+1) = node.Re;
        Pr_air_seg(k+1) = node.Pr;  K_seg(k+1)      = node.K;
        f_air_seg(k+1)  = node.f;
        T_cool_seg(k+1) = T_cool(k+1);

        % Update shared scalar outputs
        v_channel_coolant = cool_k.v;
        Re_coolant        = cool_k.Re;
        h_coolant         = cool_k.h;
        h_air             = h_air_k;
        NTU               = NTU_k;
        R_tot             = R_tot_k;
        UA_unit           = UA_k;

    end

    T_cool_seg(N_segments+1) = T_cool(N_segments+1);
    T_c_o_disc = T_air_seg(end);
    dp_hx_disc = dp_total;

    fprintf("Q_total = %.2f W,  T_air_out = %.2f K\n", Q_total, T_air_seg(end));
    % fprintf("Coolant outlet T_cool(1) = %.2f K,  T_h_o requirement = %.2f K,  deviation = %.4f K\n", ...
    %         T_cool(1), T_h_o, T_cool(1) - T_h_o);
    % fprintf("T_cool distribution: %f", T_cool_seg);
    fprintf("T_cool distribution (K): ");
    fprintf("%.4f  ", T_cool_seg);
    fprintf("\n");

end

%% ---- Length iteration ---- %%

fprintf("Total heat loss per module = %f W\n", Q_tot/n_modules);

L_low  = 1e-3;
L_high = 1;

while true
    Q_high = Q_pred_L_disc(L_high);
    fprintf("L bracket upper = %.4f m,  Q = %.2f W\n", L_high, Q_high);
    if Q_high >= (Q_tot/n_modules)
        break
    end
    L_high = L_high*2;
    if L_high > 20
        fprintf("Max heat loss possible at L = %.4f m is Q = %.2f W\n", ...
                L_high/2, Q_pred_L_disc(L_high/2));
        error('L_high grew too large — check inputs or UA_unit');
    end
end

L_solution = fzero(@(L) Q_pred_L_disc(L) - (Q_tot/n_modules), [L_low, L_high]);
fprintf("Effective length of HX = %.4f m\n", L_solution);

Q_pred_L_disc(L_solution);

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

v4      = M_dot_3/(pi*d3*d3*rho_air(P4,T4)/4);
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
mom_term  = M_dot_3*(-v3+v4);
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

counter = counter+1;
if counter > 200
    error('exiting HX');
end

N_fin_coolant   = 2*L_solution/(2*t_fin/sqrt(2)+b_t_coolant);
f_coolant       = (1/(-1.8*log10((0.0015/d_h_coolant)^1.11 + ...
                  (6.9/Re_coolant))))^2;
dp_coolant_loop = f_coolant*b_hx*(0.5*1082*v_channel_coolant^2)/d_h_coolant;
fprintf("Coolant side pressure loss = %.4f bar\n", dp_coolant_loop/1e5);

end