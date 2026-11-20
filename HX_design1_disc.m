function [dp_coolant_loop, d_h_air, M_dot_4, b_t_air, b_t_coolant, dp_hx, ...
          N_fin_air, N_fin_coolant, N_air_pass, N_coolant_pass, NTU, R_tot, ...
          v_channel_air, v_cool_seg, d_h_coolant, A_o_coolant, A_o_air, ...
          Re_air, Re_cool_seg, h_air, h_coolant, L_solution, UA_unit, v4, P4_0, ...
          M4, T4, T4_0, P4, F_drag_hx, M_hx, A4, drag_HX, T_cool_seg, dp_cool_seg, ...
          T_air_seg, P_air_seg, Re_air_seg, Pr_air_seg, v_air_seg, ...
          K_seg, Cf_air_seg, Nu_seg_arr, h_air_seg, eta_fin_seg, ...
          UA_seg_arr, NTU_seg_arr, eps_seg_arr, Q_seg_arr, dp_seg_arr, ...
          Q_pred_solution, T_h_o_solution, P_0_air_seg, T_0_air_seg, M_air_seg, f_air_seg, inlet_dp, outlet_dp, ...
          delta_BL_seg, A_free_seg, d_h_bulk_seg, T_mean_c_arr, T_mean_h_arr, ...
          h_cool_seg, Pr_cool_seg, Cf_coolant_seg, Nu_cool_seg, f_coolant_seg] = ...
          HX_design1_disc(use_DNS, e, r, hx_theta, counter, A3, v3, R, P3, ...
          d3, T_h_o, n_modules, T_h_i, T_c_i, T_mean_h, ...
          Q_tot, M_dot_coolant, M_dot_3, T3, N_segments, N_cool_seg, tol_T)
% HX_design1_disc  Discretised heat-exchanger sizing model.
%
% N_cool_seg selects the coolant-side treatment:
%   N_cool_seg == 0  ->  1D. One coolant "lane" per air segment (parallel
%                        manifold: coolant arrives at T_h_i at every air
%                        segment). Each segment iterates (inner loop) on the
%                        mean coolant temperature until the coolant outlet
%                        temperature stops changing. Coolant properties and
%                        the length-averaged Nu are evaluated at that mean.
%   N_cool_seg >= 1  ->  2D. The coolant side is split into N_cool_seg
%                        columns and marched directly (no inner loop), with
%                        local coolant properties and the LOCAL Nu at the end
%                        of each column. N_cool_seg = 1 is not the same as 0.
% Both use the same code in Q_pred_L: the cell loop runs over l = 1:N_seg_cool
% (a single cell for 1D) and the inner loop runs once for 2D.
%
% All outputs are array-shaped: (N_segments+1) x N_seg_cool for node arrays,
% N_segments x N_seg_cool for segment arrays (N_seg_cool = 1 for 1D).
%
% Friction is evaluated only inside the cell evaluation, at the cell's inlet
% state, for air (correlation or DNS) and coolant alike: Cf_air_seg, f_air_seg,
% Cf_coolant_seg, f_coolant_seg are N_segments x N_seg_cool segment arrays.
% Friction naming (outputs): Cf_air_seg / Cf_coolant_seg = Fanning (skin
% friction coefficient); f_air_seg / f_coolant_seg = Darcy friction factor
% (the one used for the pressure drops). v_air_seg = air velocity in the
% channels. Nu/friction correlations are the small helper functions below;
% the four "selection" functions (air_friction, air_Nu, coolant_friction,
% coolant_Nu) are the only place where HX_design1_disc differs from its sibling file.
%
% Nu_gn takes local = true (local Nu at position x, factor 1/3) or false
% (length-averaged over x). The air side always uses local = true at x = dx*k;
% the coolant uses local = true in 2D at x = dy*l and false in 1D at x = b_hx.
% Segment order of the air side keeps the old convention: x = dx*k.

use_2d     = N_cool_seg > 0;
N_seg_cool = max(N_cool_seg, 1);   % column count (1 for 1D)

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
                sqrt(((b_t_coolant*0.5)+(t_fin/sqrt(2)))^2 + (ht_coolant+(t_fin/sqrt(2)))^2);

A_p_air      = (b_t_air*0.5) + ((b_t_air*0.5)+(t_fin/sqrt(2)));
A_f_air      = sqrt(ht_air^2+(b_t_air*0.5)^2) + ...
               sqrt(((b_t_air*0.5)+(t_fin/sqrt(2)))^2 + (ht_air+(t_fin/sqrt(2)))^2);

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


%% ---- Shared variables written by nested functions ---- %%

T_c_o_disc        = T3;
dp_hx_disc        = 0;
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
v_air_seg    = zeros(N_segments+1, N_seg_cool);   % velocity in the channels
rho_air_seg  = zeros(N_segments+1, N_seg_cool);
M_air_seg    = zeros(N_segments+1, N_seg_cool);

T_cool_seg = zeros(N_segments, N_seg_cool+1);

% Segment arrays: N_segments x N_seg_cool
K_seg          = zeros(N_segments, N_seg_cool);
Cf_air_seg     = zeros(N_segments, N_seg_cool);   % Fanning, used in Nu (at the cell's inlet node)
f_air_seg      = zeros(N_segments, N_seg_cool);   % Darcy, used for dp_seg (at the cell's inlet node)
Q_seg_arr      = zeros(N_segments, N_seg_cool);
dp_seg_arr     = zeros(N_segments, N_seg_cool);
Nu_seg_arr     = zeros(N_segments, N_seg_cool);
h_air_seg      = zeros(N_segments, N_seg_cool);
eta_fin_seg    = zeros(N_segments, N_seg_cool);
UA_seg_arr     = zeros(N_segments, N_seg_cool);
NTU_seg_arr    = zeros(N_segments, N_seg_cool);
eps_seg_arr    = zeros(N_segments, N_seg_cool);
v_cool_seg     = zeros(N_segments, N_seg_cool);
dp_cool_seg    = zeros(N_segments, N_seg_cool);
Re_cool_seg    = zeros(N_segments, N_seg_cool);
h_cool_seg     = zeros(N_segments, N_seg_cool);
Pr_cool_seg    = zeros(N_segments, N_seg_cool);
Cf_coolant_seg = zeros(N_segments, N_seg_cool);   % Fanning
f_coolant_seg  = zeros(N_segments, N_seg_cool);   % Darcy (used for dp_cool)
Nu_cool_seg    = zeros(N_segments, N_seg_cool);

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

%% ---- Correlation helpers (physical inputs only) ---- %%
% Re, Pr: Reynolds and Prandtl number; d_h: hydraulic diameter; x: distance from the
% start of the heated length; K: temperature ratio factor (K^0.45); f_D: Darcy factor.

function f = f_lam(Re)                 % Darcy, Konakov et al. as cited in Mortean et al. 2019
    f = ((1.8*log10(Re))-1.5)^(-2);
end

function f = f_turb(Re, sr, d_h)           % Darcy, Gerl et al. 2025
    f = (1/(-2*log10(2.7*log10(Re)^1.2/Re + (sr/d_h)/3.71)))^2;
end

function f = f_haaland(Re, sr, d_h)        % Darcy, the coolant formula of the original disc model
    f = (1/(-1.8*log10((sr/d_h)^1.11 + (6.9/Re))))^2;
end

function Nu = Nu_lam(Re, Pr, d_h, x)   % laminar, length dependent (Stephan and Preusser)
    Nu = 4.364 + ((0.086)*(Re*Pr*d_h/x)^(1.33)) / (1+(0.1*Pr*(d_h*Re/x)^(0.83)));
end

function Nu = Nu_gn(Re, Pr, f_D, d_h, x, K, local)   % turbulent, Gnielinski
    if local
        c = 1/3;       % local Nu at position x
    else
        c = 1;         % Nu averaged over the length x
    end
    Nu = (f_D/8)*(Re-1000)*Pr*(1+c*(d_h/x)^(2/3))*K / ...
         (1+12.7*sqrt(f_D/8)*(Pr^(2/3)-1));
end

function Nu = Nu_DB(Re, Pr)            % turbulent, Dittus-Boelter (original disc coolant)
    Nu = 0.023*Re^0.8*Pr^0.4;
end

%% ---- Correlation selection (the only part that differs from HX_design1_disc_test) ---- %%
% BASE: hard switches, as in the original model. Air: friction used in Nu is Konakov (Fanning),
% friction used for dp is Gerl (Darcy); Nu = 3.66 (laminar), linear blend 2300-3000, Gnielinski.
% Coolant: Haaland-type friction; Nu = 4.36 (laminar), Dittus-Boelter above Re 2300.
% Cf: Fanning (used in Nu and reported); f_D: Darcy (used for pressure drops).

function [Cf, f_D] = air_friction(Re, sr, d_h)
    Cf  = f_lam(Re)/4;
    f_D = f_turb(Re, sr, d_h);
end

function f_D = coolant_friction(Re, sr, d_h)
    f_D = f_haaland(Re, sr, d_h);
end

function Nu = air_Nu(Re, Pr, f_D, d_h, x, K, local)   % f_D: Darcy (Konakov in this file)
    if Re < 2300
        Nu = 3.66;
    elseif Re < 3000
        wn = (Re-2300)/(3000-2300);
        Nu = (1-wn)*3.66 + wn*Nu_gn(Re, Pr, f_D, d_h, x, K, local);
    else
        Nu = Nu_gn(Re, Pr, f_D, d_h, x, K, local);
    end
end

function Nu = coolant_Nu(Re, Pr, ~, ~, ~, ~, ~)
    if Re < 2300
        Nu = 4.36;
    else
        Nu = Nu_DB(Re, Pr);
    end
end

%% ---- Nested: discretised Q prediction (1D and 2D in one function) ---- %%

function [Q_total, T_h_o_predicted] = Q_pred_L(L)
    % Air flows along L, discretised into N_segments slices of dx. For every air
    % segment k and coolant column l (a single column for 1D) the heat transfer is
    % evaluated from the air inlet state of that cell:
    %   1D: inner loop on the mean coolant temperature (properties at the mean,
    %       length-averaged coolant Nu);
    %   2D: one pass with the local coolant inlet state (local coolant Nu).
    % Results are written straight into the tracking arrays.

    % One full HX evaluation at a given L = one "set" in the DNS query log.
    if use_DNS
        thermoturb_query_log('new_set');
    end

    dx = L / N_segments;
    dy = b_hx / N_seg_cool;   % == b_hx for 1D

    tol_cool = 1e-4;          % K, 1D inner loop
    max_cool = 50;
    if use_2d, n_iter = 1; else, n_iter = max_cool; end

    %% ---- Air inlet node ---- %%
    T_in0   = T3;
    P_in0   = P3 - inlet_dp;
    rho_in0 = rho_air(P_in0, T_in0);
    v_in0   = (m_dot_pass_air/N_fin_air) / (rho_in0*A_o_air);
    Re_in0  = rho_in0*v_in0*d_h_air / mu_air(T_in0);
    Pr_in0  = mu_air(T_in0)*cp_air(T_in0)/k_air(T_in0);
    M_in0   = v_in0/sqrt(gamma_air(cp_air(T_in0))*R*T_in0);

    T_air_seg(1,:)  = T_in0;     P_air_seg(1,:)   = P_in0;
    v_air_seg(1,:)  = v_in0;     Re_air_seg(1,:)  = Re_in0;
    Pr_air_seg(1,:) = Pr_in0;    rho_air_seg(1,:) = rho_in0;
    M_air_seg(1,:)  = M_in0;
    P_0_air_seg(1,:) = P_in0*((1+((gamma_air(cp_air(T_in0))-1).*M_in0.^2/2)).^ ...
        ((gamma_air(cp_air(T_in0)))./(gamma_air(cp_air(T_in0))-1)));
    T_0_air_seg(1,:) = T_in0*(1+((gamma_air(cp_air(T_in0))-1)*M_in0.^2)/2);

    T_cool_seg(:,1) = T_h_i;
    T_cool_out = T3;          % 1D warm start for the coolant outlet temperature

    %% ---- Geometry that does not depend on the cell ---- %%
    N_fin_cool    = 2*L/(2*t_fin/sqrt(2)+b_t_coolant);
    A_ht_air_seg  = (((A_f_air*N_fin_air*N_air_pass) + (A_p_air*N_fin_air*N_air_pass))*dx/N_seg_cool);
    A_ht_cool_seg = ((A_f_coolant*N_fin_cool*N_coolant_pass) + ...
                     (A_p_coolant*N_fin_cool*N_coolant_pass)) * dy / N_segments;

    %% ---- Segment loop ---- %%
    for k = 1:N_segments
        if use_DNS
            thermoturb_query_log('new_seg');
        end

        for l = 1:N_seg_cool

            % Air state entering this cell (column l of the node arrays)
            T_in = T_air_seg(k,l);     P_in = P_air_seg(k,l);
            rho_in = rho_air_seg(k,l); v_in = v_air_seg(k,l);
            Re_in = Re_air_seg(k,l);   Pr_in = Pr_air_seg(k,l);

            if use_2d
                T_cool_in = T_cool_seg(k,l);   % local coolant inlet of column l
            else
                T_cool_in = T_h_i;             % parallel manifold
            end
            x_air  = dx*k;                      % distance from the air inlet
            x_cool = dy*l;                      % coolant path length to the end of column l (b_hx for 1D)
            x_air_mid = dx*(k-0.5);
            x_cool_mid = dy*(l-0.5);

            for inner = 1:n_iter

                if use_DNS
                    if use_2d
                        thermoturb_query_log('new_col');
                    else
                        thermoturb_query_log('new_inner');
                    end
                end

                if use_2d
                    T_cool_ref = T_cool_in;
                else
                    T_cool_ref = (T_h_i + T_cool_out) / 2;   % mean coolant temperature
                end

                %% -- coolant side -- %%
                rho_cool = rho_EG(T_cool_ref);
                v_cool   = (m_dot_pass_cool/N_fin_cool) / (rho_cool*A_o_coolant);
                Re_cool  = rho_cool*v_cool*d_h_coolant / mu_EG(T_cool_ref);
                Pr_cool  = mu_EG(T_cool_ref)*cp_EG_50_50(T_cool_ref)/k_EG(T_cool_ref);
                % Pr_cool_wall = 
                f_cool   = coolant_friction(Re_cool, sr, d_h_coolant);   % Darcy
                dp_cool  = f_cool*dy*(0.5*rho_cool*v_cool^2)/d_h_coolant;
                K_cool = 1; %(Pr_cool/Pr_cool_wall)^0.11
                Nu_cool  = coolant_Nu(Re_cool, Pr_cool, f_cool, d_h_coolant, x_cool_mid, K_cool, use_2d);
                h_cool   = Nu_cool*k_EG(T_cool_ref)/d_h_coolant;
                m_cool       = sqrt(2*h_cool/(k_fin_val*t_fin));
                eta_fin_cool = tanh(m_cool*ht_coolant)/(m_cool*ht_coolant);
                R_cool_seg   = 1/(h_cool*A_ht_cool_seg*eta_fin_cool);

                %% -- air side -- %%
                K_ratio = (T_in/T_cool_ref)^0.45;
                C_c     = m_dot_air/N_seg_cool*cp_air(T_in);
                C_h     = M_dot_coolant/N_segments*cp_EG_50_50(T_cool_ref);
                C_min   = min(C_c, C_h);
                C_star  = C_min/max(C_c, C_h);

                % Boundary layer analysis
                delta_BL = 0.37*x_air/(rho_in*v_in*x_air/mu_air(T_in))^(0.2);
                b_inner  = max(0, b_t_air -2*sqrt(3)*delta_BL);
                A_free   = sqrt(3)/4 * b_inner^2;
                d_h_bulk = d_h_air;

                % Friction and Nu are both evaluated here, at this cell's inlet node:
                % Cf_in = Fanning (used in Nu), f_in = Darcy (used for dp_seg).
                if use_DNS
                    [~, Nu_air, Cf_in] = thermoturb_cached(Re_in, Pr_in, T_in, T_cool_ref, d_h_bulk, x_air_mid);
                    f_in = 4*Cf_in;
                else
                    [Cf_in, f_in] = air_friction(Re_in, sr, d_h_bulk);
                    Nu_air = air_Nu(Re_in, Pr_in, 4*Cf_in, d_h_bulk, x_air_mid, K_ratio, true);   % Nu takes Darcy
                end
                h_air       = Nu_air*k_air(T_in)/d_h_bulk;
                m_air       = sqrt(2*h_air/(k_fin_val*t_fin));
                eta_fin_air = tanh(m_air*ht_air)/(m_air*ht_air);
                R_air       = 1/(h_air*A_ht_air_seg*eta_fin_air);

                %% -- heat transfer -- %%
                R_tot   = R_air + R_cond_seg + R_cool_seg;
                UA_unit = 1/R_tot;
                NTU     = UA_unit/C_min;
                eps     = epsilon_counterflow(NTU, C_star);

                Q_seg  = eps*C_min*(T_cool_in - T_in);
                dp_seg = f_in*(dx/d_h_air)*0.5*rho_in*v_in^2;

                % Outlet node conditions
                T_out = T_in + Q_seg/C_c;
                P_out = P_in - dp_seg;
                T_cool_out_new = T_cool_in - Q_seg/C_h;

                rho_out = rho_air(P_out, T_out);
                v_out   = (m_dot_pass_air/N_fin_air) / (rho_out*A_o_air);
                Re_out  = rho_out*v_out*d_h_bulk / mu_air(T_out);
                Pr_out  = mu_air(T_out)*cp_air(T_out)/k_air(T_out);

                if use_2d
                    break                                   % direct march: no inner iteration
                end
                difference = T_cool_out_new - T_cool_out;
                T_cool_out = T_cool_out_new;
                if abs(difference) < tol_cool
                    break
                elseif inner == max_cool
                    fprintf("Segment %.0f did not converge, T_cool_out: %.2f K, diff: %.5f  ", k, T_cool_out, difference)
                end
            end

            %% -- store this cell -- %%
            if use_2d
                T_mean_h_seg = (T_cool_in + T_cool_out_new) / 2;
            else
                T_mean_h_seg = T_cool_ref;
            end
            if T_out == T_in || T_mean_h_seg == T_out || T_mean_h_seg == T_in
                T_mean_c_seg = T_mean_h_seg;               % degenerate LMTD (near-zero Q_seg)
            else
                T_lm_seg = ((T_mean_h_seg-T_out)-(T_mean_h_seg-T_in)) / ...
                           log((T_mean_h_seg-T_out)/(T_mean_h_seg-T_in));
                T_mean_c_seg = T_mean_h_seg - T_lm_seg;
            end
            T_mean_c_arr(k,l) = T_mean_c_seg;
            T_mean_h_arr(k,l) = T_mean_h_seg;

            Q_seg_arr(k,l)   = Q_seg;      dp_seg_arr(k,l)  = dp_seg;
            Nu_seg_arr(k,l)  = Nu_air;     h_air_seg(k,l)   = h_air;
            eta_fin_seg(k,l) = eta_fin_air;
            UA_seg_arr(k,l)  = UA_unit;    eps_seg_arr(k,l) = eps;
            NTU_seg_arr(k,l) = NTU;        K_seg(k,l)       = K_ratio;

            v_cool_seg(k,l)  = v_cool;     dp_cool_seg(k,l)    = dp_cool;
            Re_cool_seg(k,l) = Re_cool;    h_cool_seg(k,l)     = h_cool;
            Pr_cool_seg(k,l) = Pr_cool;    Nu_cool_seg(k,l)    = Nu_cool;
            f_coolant_seg(k,l)  = f_cool;  Cf_coolant_seg(k,l) = f_cool/4;

            delta_BL_seg(k,l) = delta_BL;  A_free_seg(k,l) = A_free;
            d_h_bulk_seg(k,l) = d_h_bulk;

            T_cool_seg(k,l+1) = T_cool_out_new;

            % Air outlet node = inlet of segment k+1
            T_air_seg(k+1,l)  = T_out;     P_air_seg(k+1,l)   = P_out;
            v_air_seg(k+1,l)  = v_out;     Re_air_seg(k+1,l)  = Re_out;
            Pr_air_seg(k+1,l) = Pr_out;    rho_air_seg(k+1,l) = rho_out;
            Cf_air_seg(k,l) = Cf_in;       f_air_seg(k,l)   = f_in;    % evaluated in this cell

            M_out = v_out/sqrt(gamma_air(cp_air(T_out))*R*T_out);
            M_air_seg(k+1,l)   = M_out;
            P_0_air_seg(k+1,l) = P_out*((1+((gamma_air(cp_air(T_out))-1).*M_out.^2/2)).^ ...
                ((gamma_air(cp_air(T_out)))./(gamma_air(cp_air(T_out))-1)));
            T_0_air_seg(k+1,l) = T_out*(1+((gamma_air(cp_air(T_out))-1)*M_out.^2)/2);

            h_coolant = h_cool;     % shared scalar output (last cell)
        end
    end

    %% ---- Accounting (reads the tracking arrays back) ---- %%

    Q_total = sum(Q_seg_arr(:));

    % mean over coolant columns (no-op for 1D), then sum over segments
    dp_hx_disc = inlet_dp + sum(mean(dp_seg_arr,2));

    T_c_o_disc = mean(T_air_seg(end,:));

    fprintf("\nQ_total = %.2f W,  T_air_out = %.2f K\n", Q_total, T_c_o_disc);

    T_h_o_predicted = mean(T_cool_seg(:,end));

    fprintf("Predicted mixed coolant outlet = %.4f K,  T_h_o requirement = %.4f K,  deviation = %.4f K\n", ...
        T_h_o_predicted, T_h_o, T_h_o_predicted - T_h_o);

end

function Q = Q_only(L)
    [Q, ~] = Q_pred_L(L);
end

function T = T_only(L)
    [~, T] = Q_pred_L(L);
end

%% ---- Length iteration (shared) ---- %%

fprintf("Total heat loss per module = %f W\n", Q_tot/n_modules);

L_low  = 1e-3;
L_high = 0.5;

while true

    [Q_high, T_high] = Q_pred_L(L_high);
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
        [Q_best, T_best] = Q_pred_L(L_high/2);
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
[Q_pred_solution, T_h_o_solution] = Q_pred_L(L_solution);

%% ---- Estimate Outlet pressure loss (Kays & London, 1960) ---- %%
% Array-mean form; for N_cool_seg == 0 (N_seg_cool == 1) numerically
% identical to indexing the single column directly.

Ke = kays_london_triangular(sigma, mean(Re_air_seg(end,:)), 'Ke', false);
outlet_dp = mean(v_air_seg(end,:))^2 * rho_air(mean(P_air_seg(end,:)), mean(T_air_seg(end,:))) / 2 * (1 - sigma^2 - Ke);

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
                 sqrt(((b_t_air*0.5)+(t_fin/sqrt(2)))^2 + (ht_air+(t_fin/sqrt(2)))^2);
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
% sum(dp_cool_seg,2) then mean collapses to mean(dp_cool_seg) when
% N_seg_cool == 1 -- no if/else needed.
dp_cool_totals  = sum(dp_cool_seg,2);
dp_coolant_loop = mean(dp_cool_totals);
fprintf("Coolant side pressure loss = %.4f bar\n", dp_coolant_loop/1e5);

end