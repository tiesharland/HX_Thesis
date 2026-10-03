function [dp_coolant_loop, d_h_air, M_dot_4, b_t_air, b_t_coolant, dp_hx, ...
          N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, NTU, R_tot, ...
          v_channel_air, v_cool_seg, d_h_coolant, A_o_coolant, A_o_air, ...
          Re_air, Re_cool_seg, h_air, h_coolant, L_solution, UA_unit, v4, P4_0, ...
          M4, T4, T4_0, P4, F_drag_hx, M_hx, A4, drag_HX, T_cool_seg, dp_cool_seg, ...
          T_air_seg, P_air_seg, v_air_seg, Re_air_seg, Pr_air_seg, v_channel_seg, ...
          K_seg, f_air_seg, Nu_seg_arr, h_air_seg, eta_fin_seg, ...
          UA_seg_arr, NTU_seg_arr, eps_seg_arr, Q_seg_arr, dp_seg_arr, ...
          Q_pred_solution, T_h_o_solution, P_0_air_seg, T_0_air_seg, M_air_seg, f_hx_seg, inlet_dp, outlet_dp, ...
          delta_BL_seg, A_free_seg, d_h_bulk_seg, T_mean_c_arr, T_mean_h_arr] = ...
          HX_design1_disc(use_DNS, e, r, hx_theta, counter, A3, v3, R, P3, ...
          d3, T_h_o, n_modules, T_h_i, T_c_i, T_mean_h, ...
          Q_tot, M_dot_coolant, M_dot_3, T3, N_segments, N_cool_seg, tol_T)
% HX_design1_disc  Discretised heat-exchanger sizing model.
%
% N_cool_seg dispatches between two genuinely different algorithms for the
% coolant side (parallel-manifold coolant architecture: coolant arrives
% fresh at T_h_i at every air-flow segment position, which is why both
% paths can share the same array-shaped output layout below):
%
%   N_cool_seg == 0  ->  OLD ITERATIVE 1D algorithm. One coolant "lane"
%                        per air segment; each segment runs its own inner
%                        convergence loop on the segment's mean coolant
%                        temperature (T_mean_h_seg) until T_cool_out
%                        stabilises.
%
%   N_cool_seg >= 1  ->  2D DIRECT-MARCHING algorithm. The coolant side is
%                        split into N_cool_seg width-wise columns; the
%                        solution marches directly segment-by-segment and
%                        column-by-column with NO inner convergence loop.
%                        N_cool_seg = 1 is NOT equivalent to N_cool_seg ==
%                        0 -- it is the 2D algorithm trivially run with a
%                        single coolant column, and behaves differently
%                        from the 1D path (e.g. it uses the column's
%                        local inlet coolant temperature T_cool_seg(k,l)
%                        rather than an iterated segment-mean temperature).
%
% All outputs are returned in the array-shaped (2D-style) layout -
% (N_segments+1) x N_seg_cool for node arrays, N_segments x N_seg_cool for
% segment arrays - with N_seg_cool = 1 when N_cool_seg == 0, since the 1D
% path naturally produces single-column arrays under the parallel-manifold
% assumption.
%
% N_seg_cool == 1 collapse: every "/N_seg_cool" formula below (area
% scaling, flow length, C_c's air-mass-flow split) reduces EXACTLY to the
% 1D path's original, unscaled formula when N_seg_cool == 1. That is used
% throughout to share formulas between paths without an if/else: dy ==
% b_hx, 1/N_seg_cool == 1, m_dot_air/N_seg_cool == m_dot_air, etc. Only
% where the two paths are genuinely different ALGORITHMS (the inner
% coolant convergence loop vs. direct marching over columns) do they stay
% as two separate functions.
%
% Per-segment physics is shared between the two paths through two local
% functions, calculate_air_side and calculate_coolant_side, so a change to
% either correlation only needs to be made once:
%   - calculate_air_side(...)     Nu, h, fin efficiency and resistance on
%                                  the air side from local inlet conditions.
%   - calculate_coolant_side(...) coolant velocity/Re/Pr, h, dp, fin
%                                  efficiency, area and resistance from a
%                                  coolant reference temperature, a flow
%                                  length and an area scale factor.
% calculate_segment_1d / calculate_segment_2d do one segment's (or one
% (k,l) cell's) full physics using those two functions and write directly
% into the shared tracking arrays. Q_pred_L_disc itself just: picks the
% initial state, runs a single `for k = 1:N_segments` loop dispatching to
% whichever calculate_segment_* applies, then does the accounting
% (Q_total, dp_hx_disc, T_c_o_disc, T_h_o_predicted) by reading those same
% arrays back -- no separate 1D/2D wrapper functions.
%
% Mean-temperature outputs (T_mean_c_arr/T_mean_h_arr): both paths
% populate these the same way -- a coolant-segment mean temperature
% (1D: the converged T_mean_h_seg from its inner iteration; 2D: the
% arithmetic mean of the column's local inlet/outlet coolant temperature,
% (T_in_cool+T_out_cool)/2) fed through the same LMTD-based air-side mean
% calculation. Previously 2D hardcoded these to 0.
%
% Known, deferred issue in the N_cool_seg >= 1 (2D) path only (left as in
% HX_design1_disc_2d_v2.m; not changed by this refactor):
%   - A_ht_cool_seg's area scaling (1/N_seg_cool) and dp_cool's flow
%     length (dy) are each path's own existing choice, carried through
%     calculate_coolant_side's area_scale/flow_length arguments unchanged
%     -- not reconciled or fixed here.

use_2d     = N_cool_seg > 0;
N_seg_cool = max(N_cool_seg, 1);   % column count used for array sizing
                                    % AND for the "/N_seg_cool" collapse
                                    % described above (== 1 for 1D).

%% ---- Geometry and passage counts (shared by both algorithms) ---- %%

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
m_dot_pass_cool    = M_dot_coolant / N_coolant_pass;

A_p_coolant   = (b_t_coolant*0.5) + ((b_t_coolant*0.5)+(t_fin/sqrt(2)));
A_f_coolant   = sqrt(ht_coolant^2+(b_t_coolant*0.5)^2) + ...
    sqrt(((b_t_coolant*0.5)+(t_fin/sqrt(2)))^2 + ...
    (ht_coolant+(t_fin/sqrt(2)))^2);

A_p_air      = (b_t_air*0.5) + ((b_t_air*0.5)+(t_fin/sqrt(2)));
A_f_air      = sqrt(ht_air^2+(b_t_air*0.5)^2) + ...
    sqrt(((b_t_air*0.5)+(t_fin/sqrt(2)))^2 + ...
    (ht_air+(t_fin/sqrt(2)))^2);

% Collapses to the 1D path's original (unscaled) formula when
% N_seg_cool == 1 -- no if/else needed.
R_cond_seg = t_plate/(k_fin_val*t_plate*(2*N_pass_tot+2)*b_hx/N_seg_cool);

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

sigma = A_o_air*N_fin_air*N_air_pass/A_frontal;
fprintf('Sigma: %f\n', sigma)
Kc = kays_london_triangular(sigma, Re_air, 'Kc', false);
inlet_dp = v_channel_air^2 * rho_air(P3, T3) / 2 * (1 - sigma^2 + Kc);

fprintf('Entrance pressure drop: %.2f Pa\n', inlet_dp)

%% ---- Shared variables written by nested function ---- %%

T_c_o_disc        = T3;
dp_hx_disc        = 0;
Re_coolant        = 0;
h_air             = 0;
h_coolant         = 0;
NTU               = 0;
R_tot             = 0;
UA_unit           = 0;

% Node arrays: (N_segments+1) x N_seg_cool
T_air_seg    = zeros(N_segments+1, N_seg_cool);
T_0_air_seg  = zeros(N_segments+1, N_seg_cool);
P_air_seg    = zeros(N_segments+1, N_seg_cool);
P_0_air_seg  = zeros(N_segments+1, N_seg_cool);
Re_air_seg   = zeros(N_segments+1, N_seg_cool);
Pr_air_seg   = zeros(N_segments+1, N_seg_cool);
f_air_seg    = zeros(N_segments+1, N_seg_cool);
v_channel_seg    = zeros(N_segments+1, N_seg_cool);
v_air_seg    = zeros(N_segments+1, N_seg_cool);
rho_air_seg  = zeros(N_segments+1, N_seg_cool);
M_air_seg    = zeros(N_segments+1, N_seg_cool);
f_hx_seg     = zeros(N_segments+1, N_seg_cool);

T_cool_seg = zeros(N_segments, N_seg_cool+1);

% Segment arrays: N_segments x N_seg_cool
K_seg        = zeros(N_segments, N_seg_cool);
Q_seg_arr    = zeros(N_segments, N_seg_cool);
dp_seg_arr   = zeros(N_segments, N_seg_cool);
Nu_seg_arr   = zeros(N_segments, N_seg_cool);
h_air_seg    = zeros(N_segments, N_seg_cool);
eta_fin_seg  = zeros(N_segments, N_seg_cool);
UA_seg_arr   = zeros(N_segments, N_seg_cool);
NTU_seg_arr  = zeros(N_segments, N_seg_cool);
eps_seg_arr  = zeros(N_segments, N_seg_cool);
v_cool_seg   = zeros(N_segments, N_seg_cool);
dp_cool_seg  = zeros(N_segments, N_seg_cool);
Re_cool_seg  = zeros(N_segments, N_seg_cool);
h_cool_seg   = zeros(N_segments, N_seg_cool);

delta_BL_seg     = zeros(N_segments, N_seg_cool);
A_free_seg       = zeros(N_segments, N_seg_cool);
d_h_bulk_seg     = zeros(N_segments, N_seg_cool);
T_mean_c_arr     = zeros(N_segments, N_seg_cool);
T_mean_h_arr     = zeros(N_segments, N_seg_cool);

%% ---- Nested: counter-flow effectiveness ---- %%

function eps = epsilon_counterflow(NTU, C_star)
    if C_star < 1e-3
        eps = 1 - exp(-NTU);
    else
        eps = (1 - exp(-NTU.*(1-C_star))) ./ ...
              (1 - C_star.*exp(-NTU.*(1-C_star)));
    end
end

%% ---- Nested: shared air-side sizing (used by both 1D and 2D) ---- %%
% T_in/Re_in/Pr_in/f_in: local air inlet state. K_ratio: the T_hot/T_in
% ratio used by the Nu correlation (T_mean_h_seg/T_in for 1D,
% T_in_cool/T_in for 2D). k_idx: the CURRENT k loop index -- note this
% reproduces a pre-existing quirk shared by both original source files,
% where the Nu correlation's "(dx*k)" term actually uses the segment
% INDEX k, not a thermal conductivity, despite the name; preserved as-is,
% not fixed here. T_hot_ref: reference hot-side temperature passed to the
% DNS lookup (T_mean_h_seg for 1D, T_in_cool for 2D).

function [Nu_seg, h_air_loc, eta_fin_loc, R_air_loc, f_out_dns] = ...
        calculate_air_side(T_in, Re_in, Pr_in, f_in, K_ratio, k_idx, ...
                            d_h_bulk, dx, A_ht_air_seg, T_hot_ref)
    f_out_dns = [];
    if use_DNS
        % DNS lookup: Nu and Cf from thermoturb table
        [~, Nu_seg, Cf_seg] = thermoturb_cached(Re_in, Pr_in, T_in, T_hot_ref);
        f_out_dns = 4*Cf_seg;   % Fanning -> d_h-based
    else
        % Use correlations for Nu
        if Re_in < 2300
            Nu_seg = 3.66;
        elseif Re_in < 3000
            w       = (Re_in-2300)/(3000-2300);
            Nu_turb = (f_in/2)*(Re_in-1000)*Pr_in* ...
                      (1+(d_h_bulk/(dx*k_idx))^(2/3))*K_ratio / ...
                      (1+12.7*sqrt(f_in/2)*(Pr_in^(2/3)-1));
            Nu_seg  = (1-w)*3.66 + w*Nu_turb;
        else
            Nu_seg  = (f_in/2)*(Re_in-1000)*Pr_in* ...
                      (1+(d_h_bulk/(dx*k_idx))^(2/3))*K_ratio / ...
                      (1+12.7*sqrt(f_in/2)*(Pr_in^(2/3)-1));
        end
    end

    h_air_loc   = Nu_seg*k_air(T_in)/d_h_bulk;
    m_air       = sqrt(2*h_air_loc/(k_fin_val*t_fin));
    eta_fin_loc = tanh(m_air*ht_air)/(m_air*ht_air);
    R_air_loc   = 1/(h_air_loc*A_ht_air_seg*eta_fin_loc);
end

%% ---- Nested: shared coolant-side sizing (used by both 1D and 2D) ---- %%
% T_cool_ref: coolant reference temperature properties are evaluated at
% (T_mean_h_seg for 1D, T_in_cool for 2D). N_fin_cool: fin count for this
% dx (identical formula/value in both paths, computed once per L by the
% caller). flow_length/area_scale: each path's own existing
% parameterisation of the coolant flow path length (dp_cool) and the
% heat-transfer area scaling (A_ht_cool_seg). Both callers now always pass
% dy/(1/N_seg_cool); for the 1D path (N_seg_cool == 1) that is numerically
% identical to the original b_hx/1, so no separate 1D-specific call shape
% is needed.

function [v_cool, Re_cool, Pr_cool, dp_cool, h_cool, eta_fin_cool, ...
          A_ht_cool_seg, R_cool_seg] = ...
        calculate_coolant_side(T_cool_ref, N_fin_cool, flow_length, area_scale)

    v_cool  = (m_dot_pass_cool/N_segments/(2*N_fin_cool)) / (rho_EG(T_cool_ref)*A_o_coolant);
    Re_cool = rho_EG(T_cool_ref)*v_cool*d_h_coolant / mu_EG(T_cool_ref);
    Pr_cool = mu_EG(T_cool_ref)*cp_EG_50_50(T_cool_ref)/k_EG(T_cool_ref);
    f_cool  = (1/(-1.8*log10((0.0015/d_h_coolant)^1.11 + (6.9/Re_cool))))^2;
    dp_cool = f_cool*flow_length*(0.5*1082*v_cool^2)/d_h_coolant;

    if Re_cool < 2300
        Nu_cool = 4.36;
        h_cool  = Nu_cool*k_EG(T_cool_ref)/d_h_coolant;
    else
        Nu_cool = 0.023*Re_cool^0.8*Pr_cool^0.4;
        j_cool  = Nu_cool/(Re_cool*Pr_cool^(1/3));
        G_cool  = rho_EG(T_cool_ref)*v_cool;
        h_cool  = j_cool*G_cool*cp_EG_50_50(T_cool_ref)/Pr_cool^(2/3);
    end

    m_cool       = sqrt(2*h_cool/(k_fin_val*t_fin));
    eta_fin_cool = tanh(m_cool*ht_coolant)/(m_cool*ht_coolant);

    A_ht_cool_seg = ((A_f_coolant*N_fin_cool*N_coolant_pass) + ...
        (A_p_coolant*N_fin_cool*N_coolant_pass)) * area_scale;

    R_cool_seg = 1/(h_cool*A_ht_cool_seg*eta_fin_cool);
end

%% ---- Nested: discretised Q prediction ---- %%

function [Q_total, T_h_o_predicted] = Q_pred_L_disc(L)
    % Air flows along L, discretised into N_segments equal slices of dx.
    % Coolant: 1D runs one lane per air segment with an inner convergence
    % loop (calculate_segment_1d); 2D splits the coolant side into
    % N_seg_cool columns and marches directly with no inner loop
    % (calculate_segment_2d). Node arrays track state at N+1 nodes
    % (inlet + segment outlets); segment arrays track integrated
    % quantities over each segment. Initialisation and the segment loop
    % are shared; only which calculate_segment_* runs each k differs.

    dx = L / N_segments;
    dy = b_hx / N_seg_cool;   % == b_hx when N_seg_cool == 1 (1D path)

    tol_cool = 1e-4; % K -- 1D path only
    max_cool = 50;   % 1D path only

    %% ---- Initialise air inlet node and coolant ---- %%

    T_in = T3;
    P_in = P3 - inlet_dp;
    v_air_in = v3;

    rho_in = rho_air(P_in, T_in);
    v_channel_in = (m_dot_pass_air/N_fin_air) / (rho_in*A_o_air);
    Re_in = rho_in*v_channel_in*d_h_air / mu_air(T_in);
    Pr_in = mu_air(T_in)*cp_air(T_in)/k_air(T_in);
    f_in = 0.25*((1.8*log10(Re_in))-1.5)^(-2); %Konakov et al as cited in Mortean et.al 2019
    f_hx_in = (1/(-2*log10(2.7*log10(Re_in)^1.2/Re_in+(sr/d_h_air)/3.71)))^2;
    M_in = v_channel_in/sqrt(gamma_air(cp_air(T_in))*R*T_in);

    T_air_seg(1,:)   = T_in;      P_air_seg(1,:)   = P_in;
    v_channel_seg(1,:)   = v_channel_in;      Re_air_seg(1,:)  = Re_in;
    Pr_air_seg(1,:)  = Pr_in;     f_air_seg(1,:)   = f_in;
    v_air_seg(1,:)   = v_air_in;        rho_air_seg(1,:) = rho_in;
    f_hx_seg(1,:)    = f_hx_in;

    M_air_seg(1,:) = M_in;
    P_0_air_seg(1,:) = P_in*((1+((gamma_air(cp_air(T_in))-1).*M_in.^2/2)).^ ...
        ((gamma_air(cp_air(T_in)))./(gamma_air(cp_air(T_in))-1)));
    T_0_air_seg(1,:) = T_in*(1+((gamma_air(cp_air(T_in))-1)*M_in.^2)/2);

    T_cool_seg(:,1) = T_h_i;

    T_cool_out = T3;   % 1D-only persistent warm-start; unused by 2D

    %% ---- Geometry ---- %%

    N_fin_cool   = 2*dx/(2*t_fin/sqrt(2)+b_t_coolant);
    A_ht_air_seg = 3*b_t_air*2*N_fin_air*N_air_pass*dx/N_seg_cool;  % == no /N_seg_cool when N_seg_cool==1

    %% ---- Nested: one segment's physics, OLD ITERATIVE 1D ---- %%
    % T_in/P_in/rho_in/v_channel_in/Re_in/Pr_in/f_in/f_hx_in and
    % T_cool_out are Q_pred_L_disc's own locals (captured by closure),
    % carried over from segment to segment exactly as in the original
    % single-loop version.

    function calculate_segment_1d(k)

        for inner = 1:max_cool

            T_mean_h_seg = (T_h_i + T_cool_out) / 2;

            [v_cool, Re_cool, ~, dp_cool, h_cool, ~, ~, R_cool_seg] = ...
                calculate_coolant_side(T_mean_h_seg, N_fin_cool, dy, 1/N_seg_cool);

            K_ratio = T_mean_h_seg / T_in;
            C_c     = m_dot_air/N_seg_cool*cp_air(T_in);
            C_h     = M_dot_coolant/N_segments*cp_EG_50_50(T_mean_h_seg);
            C_min   = min(C_c, C_h);
            C_star  = C_min/max(C_c, C_h);

            % Boundary layer analysis
            delta_BL = 0.37*dx*k/(rho_in*v_channel_in*dx*k/mu_air(T_in))^(0.2);
            b_inner  = max(0, b_t_air -2*sqrt(3)*delta_BL);
            A_free   = sqrt(3)/4 * b_inner^2;
            d_h_bulk = d_h_air;

            [Nu_seg, h_air, eta_fin, R_air, f_out_dns] = calculate_air_side( ...
                T_in, Re_in, Pr_in, f_in, K_ratio, k, d_h_bulk, dx, A_ht_air_seg, T_mean_h_seg);

            % Total resistance and heat transfer
            R_tot   = R_air + R_cond_seg + R_cool_seg;
            UA_unit = 1/R_tot;
            NTU     = UA_unit/C_min;
            eps     = epsilon_counterflow(NTU, C_star);

            Q_seg  = eps*C_min*(T_h_i - T_in);
            dp_seg = f_hx_in*(dx/d_h_air)*0.5*rho_in*v_channel_in^2;

            % Outlet node conditions
            T_out = T_in + Q_seg/C_c;
            P_out = P_in - dp_seg;
            T_cool_out_old = T_cool_out;
            T_cool_out = T_h_i - Q_seg/C_h;

            v_air_out = m_dot_air/(pi*d3*d3*rho_air(P_out,T_out)/4);

            rho_out = rho_air(P_out, T_out);
            v_channel_out = (m_dot_pass_air/N_fin_air) / (rho_out*A_o_air);
            Re_out = rho_out*v_channel_out*d_h_bulk / mu_air(T_out);
            Pr_out = mu_air(T_out)*cp_air(T_out)/k_air(T_out);
            if use_DNS
                f_out = f_out_dns;
            else
                f_out = 0.25*((1.8*log10(Re_out))-1.5)^(-2); %Konakov et al as cited in Mortean et.al 2019
            end
            f_hx_out = (1/(-2*log10(2.7*log10(Re_out)^1.2/Re_out+(sr/d_h_air)/3.71)))^2; % As cited in Gerl et. al 2025

            difference = T_cool_out - T_cool_out_old;
            if abs(difference) >= tol_cool
                if inner == max_cool
                    fprintf("Segment %.0f did not converge, T_cool_out: %.2f K, diff: %.5f  ", k, T_cool_out, difference)
                    break
                end
            else
                break
            end

        end

        T_lm_seg = ((T_mean_h_seg-T_out)-(T_mean_h_seg-T_in))/log((T_mean_h_seg-T_out)/(T_mean_h_seg-T_in));
        T_mean_c_seg = T_mean_h_seg-T_lm_seg;
        T_mean_c_arr(k) = T_mean_c_seg;
        T_mean_h_arr(k) = T_mean_h_seg;

        Q_seg_arr(k)   = Q_seg;
        dp_seg_arr(k)  = dp_seg;
        Nu_seg_arr(k)  = Nu_seg;
        h_air_seg(k)   = h_air;
        eta_fin_seg(k) = eta_fin;
        UA_seg_arr(k)  = UA_unit;
        eps_seg_arr(k) = eps;
        NTU_seg_arr(k) = NTU;
        K_seg(k)       = K_ratio;
        T_cool_seg(k,end) = T_cool_out;
        v_cool_seg(k)  = v_cool;
        dp_cool_seg(k)  = dp_cool;
        Re_cool_seg(k) = Re_cool;
        h_cool_seg(k)  = h_cool;
        delta_BL_seg(k) = delta_BL;
        A_free_seg(k) = A_free;
        d_h_bulk_seg(k) = d_h_bulk;

        % Compute and store outlet node at index k+1

        T_air_seg(k+1)   = T_out;      P_air_seg(k+1)   = P_out;
        v_channel_seg(k+1)   = v_channel_out;      Re_air_seg(k+1)  = Re_out;
        Pr_air_seg(k+1)  = Pr_out;     f_air_seg(k+1)   = f_out;
        v_air_seg(k+1)   = v_air_out;  rho_air_seg(k+1) = rho_out;
        f_hx_seg(k+1)    = f_hx_out;

        M_out = v_channel_out/sqrt(gamma_air(cp_air(T_out))*R*T_out);
        M_air_seg(k+1) = M_out;
        P_0_air_seg(k+1) = P_out*((1+((gamma_air(cp_air(T_out))-1).*M_out.^2/2)).^ ...
            ((gamma_air(cp_air(T_out)))./(gamma_air(cp_air(T_out))-1)));
        T_0_air_seg(k+1) = T_in*(1+((gamma_air(cp_air(T_in))-1)*M_out.^2)/2);

        % Update shared scalar outputs
        Re_coolant = Re_cool;
        h_coolant  = h_cool;

        % Carry over to next segment
        T_in   = T_out;
        P_in   = P_out;
        rho_in = rho_out;
        v_channel_in = v_channel_out;
        Re_in = Re_out;
        Pr_in = Pr_out;
        f_in  = f_out;
        f_hx_in = f_hx_out;

    end

    %% ---- Nested: one (k,l) cell's physics, 2D DIRECT MARCH ---- %%

    function calculate_segment_2d(k)

        for l = 1:N_seg_cool

            T_in = T_air_seg(k,l);      P_in = P_air_seg(k,l);
            rho_in = rho_air_seg(k,l);  v_channel_in = v_channel_seg(k,l);
            Re_in = Re_air_seg(k,l);    Pr_in = Pr_air_seg(k,l);
            f_in  = f_air_seg(k,l);     f_hx_in = f_hx_seg(k,l);

            T_in_cool = T_cool_seg(k,l);

            [v_cool, Re_cool, ~, dp_cool, h_cool, ~, ~, R_cool_seg] = ...
                calculate_coolant_side(T_in_cool, N_fin_cool, dy, 1/N_seg_cool);

            K_ratio = T_in_cool / T_in;
            C_c     = m_dot_air/N_seg_cool*cp_air(T_in);
            C_h     = M_dot_coolant/N_segments*cp_EG_50_50(T_in_cool);
            C_min   = min(C_c, C_h);
            C_star  = C_min/max(C_c, C_h);

            % Boundary layer analysis
            delta_BL = 0.37*dx*k/(rho_in*v_channel_in*dx*k/mu_air(T_in))^(0.2);
            b_inner  = max(0, b_t_air -2*sqrt(3)*delta_BL);
            A_free   = sqrt(3)/4 * b_inner^2;
            d_h_bulk = d_h_air;

            [Nu_seg, h_air, eta_fin, R_air, f_out_dns] = calculate_air_side( ...
                T_in, Re_in, Pr_in, f_in, K_ratio, k, d_h_bulk, dx, A_ht_air_seg, T_in_cool);

            % Total resistance and heat transfer
            R_tot   = R_air + R_cond_seg + R_cool_seg;
            UA_unit = 1/R_tot;
            NTU     = UA_unit/C_min;
            eps     = epsilon_counterflow(NTU, C_star);

            Q_seg  = eps*C_min*(T_in_cool - T_in);
            dp_seg = f_hx_in*(dx/d_h_air)*0.5*rho_in*v_channel_in^2;

            % Outlet node conditions
            T_out = T_in + Q_seg/C_c;
            P_out = P_in - dp_seg;
            T_out_cool = T_in_cool - Q_seg/C_h;

            v_air_out = m_dot_air/(pi*d3*d3*rho_air(P_out,T_out)/4);

            rho_out = rho_air(P_out, T_out);
            v_channel_out = (m_dot_pass_air/N_fin_air) / (rho_out*A_o_air);
            Re_out = rho_out*v_channel_out*d_h_bulk / mu_air(T_out);
            Pr_out = mu_air(T_out)*cp_air(T_out)/k_air(T_out);
            if use_DNS
                f_out = f_out_dns;
            else
                f_out = 0.25*((1.8*log10(Re_out))-1.5)^(-2); %Konakov et al as cited in Mortean et.al 2019
            end
            f_hx_out = (1/(-2*log10(2.7*log10(Re_out)^1.2/Re_out+(sr/d_h_air)/3.71)))^2; % As cited in Gerl et. al 2025

            % Mean-temperature outputs: same style of calculation 1D
            % uses, now populated here too (previously hardcoded to 0).
            T_mean_h_seg = (T_in_cool + T_out_cool) / 2;
            if T_out == T_in || T_mean_h_seg == T_out || T_mean_h_seg == T_in
                % Degenerate LMTD (near-zero Q_seg for this cell) -- guard
                % against log(1)/0 that 1D's formula doesn't need, since
                % 1D never has a near-zero Q_seg.
                T_mean_c_seg = T_mean_h_seg;
            else
                T_lm_seg = ((T_mean_h_seg-T_out)-(T_mean_h_seg-T_in)) / ...
                           log((T_mean_h_seg-T_out)/(T_mean_h_seg-T_in));
                T_mean_c_seg = T_mean_h_seg - T_lm_seg;
            end
            T_mean_c_arr(k,l) = T_mean_c_seg;
            T_mean_h_arr(k,l) = T_mean_h_seg;

            Q_seg_arr(k,l)   = Q_seg;
            dp_seg_arr(k,l)  = dp_seg;
            Nu_seg_arr(k,l)  = Nu_seg;
            h_air_seg(k,l)   = h_air;
            eta_fin_seg(k,l) = eta_fin;
            UA_seg_arr(k,l)  = UA_unit;
            eps_seg_arr(k,l) = eps;
            NTU_seg_arr(k,l) = NTU;
            K_seg(k,l)       = K_ratio;
            v_cool_seg(k,l)  = v_cool;
            dp_cool_seg(k,l)  = dp_cool;
            Re_cool_seg(k,l)  = Re_cool;
            h_cool_seg(k,l)   = h_cool;

            delta_BL_seg(k,l) = delta_BL;
            A_free_seg(k,l) = A_free;
            d_h_bulk_seg(k,l) = d_h_bulk;

            % Compute and store outlet node at index k+1

            T_cool_seg(k,l+1) = T_out_cool;

            T_air_seg(k+1,l)   = T_out;      P_air_seg(k+1,l)   = P_out;
            v_channel_seg(k+1,l)   = v_channel_out;      Re_air_seg(k+1,l)  = Re_out;
            Pr_air_seg(k+1,l)  = Pr_out;     f_air_seg(k+1,l)   = f_out;
            v_air_seg(k+1,l)   = v_air_out;  rho_air_seg(k+1,l) = rho_out;
            f_hx_seg(k+1,l)    = f_hx_out;

            M_out = v_channel_out/sqrt(gamma_air(cp_air(T_out))*R*T_out);
            M_air_seg(k+1,l) = M_out;
            P_0_air_seg(k+1,l) = P_out*((1+((gamma_air(cp_air(T_out))-1).*M_out.^2/2)).^ ...
                ((gamma_air(cp_air(T_out)))./(gamma_air(cp_air(T_out))-1)));
            T_0_air_seg(k+1,l) = T_out*(1+((gamma_air(cp_air(T_out))-1)*M_out.^2)/2);

        end

    end

    %% ---- Segment loop (shared) ---- %%

    for k = 1:N_segments
        if use_2d
            calculate_segment_2d(k);
        else
            calculate_segment_1d(k);
        end
    end

    %% ---- Accounting (shared; reads the tracking arrays back) ---- %%

    Q_total = sum(Q_seg_arr(:));

    % mean(dp_seg_arr,2): average across coolant columns for each k (a
    % no-op when N_seg_cool == 1); sum(...) then totals over k. Same
    % value either path.
    dp_hx_disc = inlet_dp + sum(mean(dp_seg_arr,2));

    T_c_o_disc = mean(T_air_seg(end,:));

    fprintf("\nQ_total = %.2f W,  T_air_out = %.2f K\n", Q_total, T_c_o_disc);

    T_h_o_predicted = mean(T_cool_seg(:,end));

    fprintf("Predicted mixed coolant outlet = %.4f K,  T_h_o requirement = %.4f K,  deviation = %.4f K\n", ...
        T_h_o_predicted, T_h_o, T_h_o_predicted - T_h_o);

end

function Q = Q_only(L)
    [Q, ~] = Q_pred_L_disc(L);
end

function T = T_only(L)
    [~, T] = Q_pred_L_disc(L);
end

%% ---- Length iteration (shared) ---- %%

fprintf("Total heat loss per module = %f W\n", Q_tot/n_modules);

L_low  = 1e-3;
L_high = 1;

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
    Q_residual = @(L) Q_only(L) - Q_tot/n_modules;
    L_solution = fzero(Q_residual, [L_low, L_high]);
    fprintf("L solved on Q residual: L = %.4f m\n", L_solution);
elseif T_close && Q_satisfied
    % Case 2a: T within tolerance of T_h_o - solve on Q residual
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

%% ---- Estimate Outlet pressure loss (Kays & London, 1960) ---- %%
% Array-mean form; for N_cool_seg == 0 (N_seg_cool == 1) numerically
% identical to indexing the single column directly.

Ke = kays_london_triangular(sigma, mean(Re_air_seg(end,:)), 'Ke', false);
outlet_dp = mean(v_channel_seg(end,:))^2 * rho_air(mean(P_air_seg(end,:)), mean(T_air_seg(end,:))) / 2 * (1 - sigma^2 - Ke);

fprintf('Exit pressure rise: %.2f Pa\n', outlet_dp)


dp_hx = dp_hx_disc - outlet_dp;
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

N_fin_coolant   = 2*L_solution/(2*t_fin/sqrt(2)+b_t_coolant);
f_coolant       = (1/(-1.8*log10((0.0015/d_h_coolant)^1.11 + (6.9/mean(Re_cool_seg,'all')))))^2; %#ok<NASGU> unused, kept for parity with both source files
% sum(dp_cool_seg,2) then mean collapses to mean(dp_cool_seg) when
% N_seg_cool == 1 -- no if/else needed.
dp_cool_totals  = sum(dp_cool_seg,2);
dp_coolant_loop = mean(dp_cool_totals);
fprintf("Coolant side pressure loss = %.4f bar\n", dp_coolant_loop/1e5);

end