function [dp_coolant_loop, d_h_air, M_dot_4, b_t_air, b_t_coolant, dp_hx, ...
          N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, NTU, R_tot, ...
          v_channel_air, v_channel_coolant, d_h_coolant, A_o_coolant, A_o_air, ...
          Re_air, Re_coolant, h_air, h_coolant, L_solution, UA_unit, v4, P4_0, ...
          M4, T4, T4_0, P4, F_drag_hx, M_hx, A4, drag_HX, T_cool_out_arr, dp_cool_seg, ...
          T_air_seg, P_air_seg, v_air_seg, Re_air_seg, Pr_air_seg, v_channel_seg, ...
          K_seg, f_air_seg, Nu_seg_arr, h_air_seg, eta_fin_seg, ...
          UA_seg_arr, NTU_seg_arr, eps_seg_arr, Q_seg_arr, dp_seg_arr, ...
          Q_pred_solution, T_h_o_solution, P_0_air_seg, T_0_air_seg, M_air_seg, f_hx_seg, inlet_dp, outlet_dp, ...
          delta_BL_seg, A_free_seg, d_h_bulk_seg, T_mean_c_arr, T_mean_h_arr, HX2D] = ...
          HX_design1_disc(use_DNS, e, r, hx_theta, counter, A3, v3, R, P3, ...
          d3, T_h_o, n_modules, T_h_i, T_c_i, T_mean_h, ...
          Q_tot, M_dot_coolant, M_dot_3, T3, N_segments, tol_T, D)

if nargin < 22
    D = 5; % default number of coolant (transverse/crossflow) segments if not supplied
    warning('HX_design1_disc:DdefaultUsed', ...
        'D (number of coolant transverse segments) not supplied - defaulting to D = 5.');
end

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

% Plates
A_plates = (N_pass_tot + 2) * t_plate * b_hx;

% Air fin solid material (fins between channels within each air passage)
% Total air passage area minus open channel area, times number of passages
A_one_air_pass = (ht_air + t_fin/sqrt(2)) * b_hx;  % approximate height of one air layer
A_fins_air_solid = (A_one_air_pass - A_o_air*N_fin_air) * N_air_pass;

% Coolant passage face area (seen end-on in frontal plane)
% The entire coolant layer height times b_hx, minus nothing since
% coolant open area is not in the air-flow direction
A_coolant_layer = (ht_coolant + t_fin/sqrt(2)) * b_hx * N_coolant_pass;

A_solid = A_plates + A_fins_air_solid + A_coolant_layer;
A_solid_alt = A_frontal - A_o_air*N_fin_air*N_air_pass;


%% ---- Inlet bulk reference (scalar outputs) ---- %%

v_channel_air = (m_dot_pass_air/N_fin_air) / (rho_air(P3,T3)*A_o_air);
Re_air        = rho_air(P3,T3)*v_channel_air*d_h_air / mu_air(T3);

%% ---- Estimate Inlet pressure loss (Kays & London, 1960) ---- %%
% - 0.5*rho_air(P4,T4)*v4^2
% inlet_dp = (0.5*rho_air(P3,T3)*v3^2)*(A_solid_alt/A_frontal);

sigma = A_o_air*N_fin_air*N_air_pass/A_frontal;
fprintf('Sigma: %f\n', sigma)
Kc = kays_london_triangular(sigma, Re_air, 'Kc', false);
inlet_dp = v_channel_air^2 * rho_air(P3, T3) / 2 * (1 - sigma^2 + Kc);

fprintf('Entrance pressure drop: %.2f Pa\n', inlet_dp)

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
T_0_air_seg  = zeros(1, N_segments+1);
P_air_seg    = zeros(1, N_segments+1);
P_0_air_seg  = zeros(1, N_segments+1);
Re_air_seg   = zeros(1, N_segments+1);
Pr_air_seg   = zeros(1, N_segments+1);
f_air_seg    = zeros(1, N_segments+1);
v_channel_seg    = zeros(1, N_segments+1);
v_air_seg    = zeros(1, N_segments+1);
rho_air_seg  = zeros(1, N_segments+1);
M_air_seg    = zeros(1, N_segments+1);
f_hx_seg     = zeros(1, N_segments+1);

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
delta_BL_seg     = zeros(1, N_segments);
A_free_seg       = zeros(1, N_segments);
d_h_bulk_seg     = zeros(1, N_segments);
T_mean_c_arr     = zeros(1, N_segments);
T_mean_h_arr     = zeros(1, N_segments);

% 2D grid arrays: rows = air segments k (1..N_segments, downstream),
% cols = coolant segments j (1..D, transverse/width direction).
% Cell (k,j) stores the OUTLET state of that cell.
T_air_grid   = zeros(N_segments, D);
T_cool_grid  = zeros(N_segments, D);
Q_grid       = zeros(N_segments, D);
h_air_grid   = zeros(N_segments, D);
h_cool_grid  = zeros(N_segments, D);
Re_cool_grid = zeros(N_segments, D);
NTU_grid     = zeros(N_segments, D);
eps_grid     = zeros(N_segments, D);
UA_grid      = zeros(N_segments, D);
dp_cool_grid = zeros(N_segments, D);
Nu_grid      = zeros(N_segments, D);
eta_fin_grid = zeros(N_segments, D);
K_grid       = zeros(N_segments, D);
v_cool_grid  = zeros(N_segments, D);

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
    dp_total = inlet_dp;

    tol_cool = 1e-4; % K
    max_cool = 50;

    %% ---- Initialise air inlet node and coolant ---- %%
    
    T_in = T3;
    P_in = P3 - inlet_dp;
    v_air_in = v3;

    T_cool_out = T3;
    rho_in = rho_air(P_in, T_in);
    % v_channel_in = v3*A3 / (A_o_air * N_fin_air * N_air_pass)
    v_channel_in = (m_dot_pass_air/N_fin_air) / (rho_in*A_o_air);
    Re_in = rho_in*v_channel_in*d_h_air / mu_air(T_in);
    Pr_in = mu_air(T_in)*cp_air(T_in)/k_air(T_in);
    f_in = 0.25*((1.8*log10(Re_in))-1.5)^(-2); %Konakov et al as cited in Mortean et.al 2019
    f_hx_in = (1/(-2*log10(2.7*log10(Re_in)^1.2/Re_in+(sr/d_h_air)/3.71)))^2;
    M_in = v_channel_in/sqrt(gamma_air(cp_air(T_in))*R*T_in);

    T_air_seg(1)   = T_in;      P_air_seg(1)   = P_in;
    v_channel_seg(1)   = v_channel_in;      Re_air_seg(1)  = Re_in;
    Pr_air_seg(1)  = Pr_in;     f_air_seg(1)   = f_in;     
    v_air_seg(1)   = v_air_in;        rho_air_seg(1) = rho_in;
    f_hx_seg(1)    = f_hx_in;

    M_air_seg(1) = M_in;
    P_0_air_seg(1) = P_in*((1+((gamma_air(cp_air(T_in))-1).*M_in.^2/2)).^ ...
        ((gamma_air(cp_air(T_in)))./(gamma_air(cp_air(T_in))-1)));
    T_0_air_seg(1) = T_in*(1+((gamma_air(cp_air(T_in))-1)*M_in.^2)/2);

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

    %% ---- Segment loop: air (k) x coolant (D) crossflow discretisation ---- %%
    %
    % Air (index k, 1..N_segments) flows downstream along L. It is split into
    % D transverse "lanes" (index j) of m_dot_air/D each; a given lane flows
    % straight through from k=1 to k=N_segments with NO lateral mixing
    % between lanes (unmixed-unmixed crossflow assumption).
    %
    % Coolant (index j, 1..D) flows across the width b_hx. For every air
    % segment k, the coolant re-enters at T_h_i (fresh from the inlet
    % manifold that runs the length of the HX) and flows straight across
    % from j=1 to j=D, independently of what happens at other k's.
    %
    % T_air_grid(k,j) / T_cool_grid(k,j) hold the OUTLET state of cell (k,j).
    % Bulk 1D arrays (T_air_seg, dp_seg_arr, etc.) are still populated - from
    % lane-averaged / summed cell values - so the rest of the function (the
    % L iteration, pressure loss, exit-state calcs) is unaffected.

    m_dot_air_lane = m_dot_air / D;

    for k = 1:N_segments

        % Air entering segment k, lane by lane (no lateral mixing)
        if k == 1
            T_air_in_row = T_in * ones(1, D);
        else
            T_air_in_row = T_air_grid(k-1, :);
        end

        % Coolant re-enters fresh at T_h_i for every air segment (manifold)
        T_cool_in_cell = T_h_i;

        Q_seg = 0;

        % Boundary layer analysis (bulk, based on segment inlet - shared by
        % all D cells since it does not depend on the transverse position)
        delta_BL = 0.37*dx*k/(rho_in*v_channel_in*dx*k/mu_air(T_in))^(0.2);
        b_inner  = max(0, b_t_air - 2*sqrt(3)*delta_BL);
        A_free   = sqrt(3)/4 * b_inner^2;
        d_h_bulk = d_h_air;

        for j = 1:D

            T_air_in_cell = T_air_in_row(j);

            T_cool_out_cell = T_cool_in_cell; % initial guess for the inner (property) iteration

            for inner = 1:max_cool

                T_mean_h_cell = (T_cool_in_cell + T_cool_out_cell) / 2;

                % ---- Coolant-side properties & resistance for this cell ----
                v_cool  = (m_dot_pass_cool_seg/(2*N_fin_cool)) / (rho_EG(T_mean_h_cell)*A_o_coolant);
                Re_cool = rho_EG(T_mean_h_cell)*v_cool*d_h_coolant / mu_EG(T_mean_h_cell);
                Pr_cool = mu_EG(T_mean_h_cell)*cp_EG_50_50(T_mean_h_cell)/k_EG(T_mean_h_cell);
                f_cool  = (1/(-1.8*log10((0.0015/d_h_coolant)^1.11 + (6.9/Re_cool))))^2;
                dp_cool_cell = f_cool*(b_hx/D)*(0.5*1082*v_cool^2)/d_h_coolant;

                if Re_cool < 2300
                    Nu_cool = 4.36;
                    h_cool  = Nu_cool*k_EG(T_mean_h_cell)/d_h_coolant;
                else
                    Nu_cool = 0.023*Re_cool^0.8*Pr_cool^0.4;
                    j_cool  = Nu_cool/(Re_cool*Pr_cool^(1/3));
                    G_cool  = rho_EG(T_mean_h_cell)*v_cool;
                    h_cool  = j_cool*G_cool*cp_EG_50_50(T_mean_h_cell)/Pr_cool^(2/3);
                end

                m_cool       = sqrt(2*h_cool/(k_fin_val*t_fin));
                eta_fin_cool = tanh(m_cool*ht_coolant)/(m_cool*ht_coolant);

                A_ht_cool_cell = ((A_f_coolant*N_fin_cool*N_coolant_pass) + ...
                                   (A_p_coolant*N_fin_cool*N_coolant_pass)) / D;

                R_cool_cell = 1/(h_cool*A_ht_cool_cell*eta_fin_cool);

                % ---- Cell capacities ----
                K_cell  = T_mean_h_cell / T_air_in_cell;
                C_c     = m_dot_air_lane*cp_air(T_air_in_cell);
                C_h     = m_dot_cool_seg*cp_EG_50_50(T_mean_h_cell);
                C_min   = min(C_c, C_h);
                C_star  = C_min/max(C_c, C_h);

                % ---- Air-side heat transfer (channel Re/Pr/f are unaffected
                % by the D-grouping - same physical channels - evaluated at
                % this cell's local air inlet temperature) ----
                if use_DNS
                    [~, Nu_cell, Cf_cell] = thermoturb_cached(Re_in, Pr_in, T_air_in_cell, T_mean_h_cell);
                    f_air_cell = 4*Cf_cell;
                else
                    f_air_cell = f_in;
                    if Re_in < 2300
                        Nu_cell = 3.66;
                    elseif Re_in < 3000
                        w        = (Re_in-2300)/(3000-2300);
                        Nu_turb  = (f_air_cell/2)*(Re_in-1000)*Pr_in* ...
                                   (1+(d_h_bulk/(dx*k))^(2/3))*K_cell / ...
                                   (1+12.7*sqrt(f_air_cell/2)*(Pr_in^(2/3)-1));
                        Nu_cell  = (1-w)*3.66 + w*Nu_turb;
                    else
                        Nu_cell  = (f_air_cell/2)*(Re_in-1000)*Pr_in* ...
                                   (1+(d_h_bulk/(dx*k))^(2/3))*K_cell / ...
                                   (1+12.7*sqrt(f_air_cell/2)*(Pr_in^(2/3)-1));
                    end
                end

                h_air_cell    = Nu_cell*k_air(T_air_in_cell)/d_h_air;
                m_air         = sqrt(2*h_air_cell/(k_fin_val*t_fin));
                eta_fin_cell  = tanh(m_air*ht_air)/(m_air*ht_air);
                A_ht_air_cell = A_ht_air_seg / D;
                R_air_cell    = 1/(h_air_cell*A_ht_air_cell*eta_fin_cell);

                % ---- Combine and solve this cell ----
                R_cond_cell = R_cond_seg * D; % same area-scaling logic as A_ht_*_cell above
                R_tot_cell  = R_air_cell + R_cond_cell + R_cool_cell;
                UA_cell     = 1/R_tot_cell;
                NTU_cell    = UA_cell/C_min;
                eps_cell    = epsilon_counterflow(NTU_cell, C_star);

                Q_cell = eps_cell*C_min*(T_cool_in_cell - T_air_in_cell);

                T_cool_out_cell_new = T_cool_in_cell - Q_cell/C_h;

                difference = T_cool_out_cell_new - T_cool_out_cell;
                T_cool_out_cell = T_cool_out_cell_new;

                if abs(difference) < tol_cool
                    break
                elseif inner == max_cool
                    fprintf("Cell (air=%.0f, cool=%.0f) did not converge, T_cool_out: %.2f K, diff: %.5f\n", ...
                        k, j, T_cool_out_cell, difference)
                end

            end

            T_air_out_cell = T_air_in_cell + Q_cell/C_c;

            % ---- Store cell results in the 2D grids ----
            T_air_grid(k,j)   = T_air_out_cell;
            T_cool_grid(k,j)  = T_cool_out_cell;
            Q_grid(k,j)       = Q_cell;
            h_air_grid(k,j)   = h_air_cell;
            h_cool_grid(k,j)  = h_cool;
            Re_cool_grid(k,j) = Re_cool;
            NTU_grid(k,j)     = NTU_cell;
            eps_grid(k,j)     = eps_cell;
            UA_grid(k,j)      = UA_cell;
            dp_cool_grid(k,j) = dp_cool_cell;
            Nu_grid(k,j)      = Nu_cell;
            eta_fin_grid(k,j) = eta_fin_cell;
            K_grid(k,j)       = K_cell;
            v_cool_grid(k,j)  = v_cool;

            Q_seg = Q_seg + Q_cell;
            T_cool_in_cell = T_cool_out_cell; % coolant advances to the next transverse sub-segment

            % Keep shared scalar outputs updated (last cell -> final value,
            % matching the original function's convention)
            v_channel_coolant = v_cool;
            Re_coolant        = Re_cool;
            h_coolant         = h_cool;
            h_air             = h_air_cell;
            NTU               = NTU_cell;
            R_tot             = R_tot_cell;
            UA_unit           = UA_cell;

        end % j (coolant / transverse) loop

        % ---- Bulk (mass-flow-weighted mean) air state for this length ----
        % position, used by the rest of the 1D length-marching model.
        T_out = mean(T_air_grid(k, :));

        dp_seg = f_hx_in*(dx/d_h_air)*0.5*rho_in*v_channel_in^2;
        P_out  = P_in - dp_seg;

        rho_out       = rho_air(P_out, T_out);
        v_channel_out = (m_dot_pass_air/N_fin_air) / (rho_out*A_o_air);
        Re_out        = rho_out*v_channel_out*d_h_bulk / mu_air(T_out);
        Pr_out        = mu_air(T_out)*cp_air(T_out)/k_air(T_out);
        if ~use_DNS
            f_out = 0.25*((1.8*log10(Re_out))-1.5)^(-2); %Konakov et al as cited in Mortean et.al 2019
        else
            f_out = f_air_cell;
        end
        f_hx_out = (1/(-2*log10(2.7*log10(Re_out)^1.2/Re_out+(sr/d_h_air)/3.71)))^2; % As cited in Gerl et. al 2025

        T_cool_out = T_cool_in_cell; % coolant leaving the LAST coolant sub-segment (j=D) for this k

        T_mean_h_seg = mean([T_h_i, T_cool_grid(k, :)]);
        T_lm_seg = ((T_mean_h_seg-T_out)-(T_mean_h_seg-T_in))/log((T_mean_h_seg-T_out)/(T_mean_h_seg-T_in));
        T_mean_c_seg = T_mean_h_seg - T_lm_seg;
        T_mean_c_arr(k) = T_mean_c_seg;
        T_mean_h_arr(k) = T_mean_h_seg;

        Q_seg_arr(k)      = Q_seg;
        dp_seg_arr(k)     = dp_seg;
        Nu_seg_arr(k)     = mean(Nu_grid(k, :));
        h_air_seg(k)      = mean(h_air_grid(k, :));
        eta_fin_seg(k)    = mean(eta_fin_grid(k, :));
        UA_seg_arr(k)     = sum(UA_grid(k, :));   % total conductance across the D cells (series path)
        eps_seg_arr(k)    = mean(eps_grid(k, :));
        NTU_seg_arr(k)    = mean(NTU_grid(k, :));
        K_seg(k)          = mean(K_grid(k, :));
        T_cool_out_arr(k) = T_cool_out;
        v_cool_seg(k)     = mean(v_cool_grid(k, :));
        dp_cool_seg(k)    = sum(dp_cool_grid(k, :)); % total coolant dp across D sub-segments in series
        delta_BL_seg(k)   = delta_BL;
        A_free_seg(k)     = A_free;
        d_h_bulk_seg(k)   = d_h_bulk;

        % Compute and store outlet node at index k+1

        T_air_seg(k+1)   = T_out;      P_air_seg(k+1)   = P_out;
        v_channel_seg(k+1)   = v_channel_out;      Re_air_seg(k+1)  = Re_out;
        Pr_air_seg(k+1)  = Pr_out;     f_air_seg(k+1)   = f_out;
        v_air_seg(k+1)   = m_dot_air/(pi*d3*d3*rho_out/4);  rho_air_seg(k+1) = rho_out;
        f_hx_seg(k+1)    = f_hx_out;

        M_out = v_channel_out/sqrt(gamma_air(cp_air(T_out))*R*T_out);
        M_air_seg(k+1) = M_out;
        P_0_air_seg(k+1) = P_out*((1+((gamma_air(cp_air(T_out))-1).*M_out.^2/2)).^ ...
            ((gamma_air(cp_air(T_out)))./(gamma_air(cp_air(T_out))-1)));
        T_0_air_seg(k+1) = T_in*(1+((gamma_air(cp_air(T_in))-1)*M_out.^2)/2);

        T_in   = T_out;
        P_in   = P_out;
        rho_in = rho_out;
        v_channel_in = v_channel_out;
        Re_in = Re_out;
        Pr_in = Pr_out;
        f_in  = f_out;
        f_hx_in = f_hx_out;

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
    dp_hx_disc = dp_total; %+ inlet_dp/2;

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
        % Case 1: both satisfied — Solve for Q %%prefer T_h_o residual
        fprintf("Case 1: both constraints satisfied at L_high = %.4f m\n", L_high);
        break

    elseif Q_satisfied && ~T_satisfied
        % Case 2: Q satisfied but T_h_o not reached
        if T_close
            fprintf('Case 2a: Q satisfied and T within tol_T at L_high = %.4f m\n', L_high);
            break
        else
            fprintf('Case 2b: Q satisfied and T not within tol_T\n');
            break
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
    % % options    = optimset('TolFun', tol_L_iter);
    % T_residual = @(L) T_only(L) - T_h_o;
    % L_solution = fzero(T_residual, [L_low, L_high]);
    % fprintf("L solved on T_h_o residual: L = %.4f m\n", L_solution);
    Q_residual = @(L) Q_only(L) - Q_tot/n_modules;
    L_solution = fzero(Q_residual, [L_low, L_high]);
    fprintf("L solved on Q residual: L = %.4f m\n", L_solution);
elseif T_close && Q_satisfied
    % Case 2a: T within tolerance of T_h_o - solve on T_h_o + tol_T residual
    % T_residual = @(L) T_only(L) - (T_h_o + tol_T);
    % L_solution = fzero(T_residual, [L_low, L_high]);
    % fprintf("L solved on T_h_o + %.1f K residual: L = %.4f m\n", tol_T, L_solution);
    % % % Case 2a: T within tolerance of T_h_o - solve on Q residual
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

%% ---- Assemble 2D (air x coolant) discretisation matrices for output ---- %%
% Rows = air segments (k = 1..N_segments, downstream/length direction)
% Cols = coolant segments (j = 1..D, transverse/width direction)
% Each HX2D.<field>(k,j) is the OUTLET state/value of cell (k,j).
HX2D.T_air      = T_air_grid;      % [K]   air outlet temperature per cell
HX2D.T_cool     = T_cool_grid;     % [K]   coolant outlet temperature per cell
HX2D.Q          = Q_grid;          % [W]   heat transferred in each cell
HX2D.h_air      = h_air_grid;      % [W/m^2K] air-side h per cell
HX2D.h_cool     = h_cool_grid;     % [W/m^2K] coolant-side h per cell
HX2D.Re_cool    = Re_cool_grid;    % [-]   coolant Reynolds number per cell
HX2D.NTU        = NTU_grid;        % [-]   NTU per cell
HX2D.eps        = eps_grid;        % [-]   effectiveness per cell
HX2D.UA         = UA_grid;         % [W/K] conductance per cell
HX2D.dp_cool    = dp_cool_grid;    % [Pa]  coolant pressure drop per cell
HX2D.N_air_seg  = N_segments;
HX2D.N_cool_seg = D;

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

%% ---- Estimate Outlet pressure loss (Kays & London, 1960) ---- %%

Ke = kays_london_triangular(sigma, Re_air_seg(end), 'Ke', false);
outlet_dp = v_channel_seg(end)^2 * rho_air(P_air_seg(end), T_air_seg(end)) / 2 * (1 - sigma^2 - Ke);

% outlet_dp = -inlet_dp/2;

fprintf('Exit pressure rise: %.2f Pa\n', outlet_dp)


dp_hx = dp_hx_disc - outlet_dp;
% dp_hx_lump = (1/(-2*log10(2.7*log10(Re_out)^1.2/Re_out+(sr/d_h_air)/3.71)))^2;
fprintf("Total pressure loss = %.2f Pa\n", dp_hx);

%% ---- Exit state variables ---- %%

A4  = A3;
T4  = T_c_o_disc;
P4  = P3 - dp_hx;
fprintf("P4 = %.2f Pa\n", P4);

if P4 < 0
    error('HX is unfeasible; Check inputs');
end

v4      = m_dot_air/(pi*d3*d3*rho_air(P4,T4)/4);
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

% fprintf('\nPlate area: %.5f m^2, fin area: %.5f m^2, together: %.5f m^2\n', A_p_air, A_f_air, A_p_air + A_f_air);
% fprintf('Solid area: %.5f m^2\n', A_solid);
% fprintf('Open area: %.5f m^2, together: %.5f m^2\n', A_o_air*N_fin_air*N_air_pass, A_solid+A_o_air*N_fin_air*N_air_pass);
% fprintf('Diffuser outlet area: %.5f m^2, HX inlet area: %.5f m^2\n', A3, A_frontal);
% fprintf('Cross sectional area of core: %.5f m^2\n', V_core/L_solution)
% % - 0.5*rho_air(P4,T4)*v4^2
% inlet_dp = (0.5*rho_air(P3,T3)*v3^2)*(A_solid_alt/A_frontal);
% fprintf('Estimate inlet bulk pressure loss: %.2f Pa\n\n', inlet_dp)

end