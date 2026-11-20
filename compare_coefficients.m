%% compare_coefficients.m
% Nusselt number and skin friction coefficient (Cf = Fanning friction factor = f_Darcy/4)
% against Re.
%
%   AIR      : the DNS table (thermoturb_cached) and any number of correlations.
%   COOLANT  : the same correlations (no DNS data exists for the coolant).
%
% The correlations are local functions at the bottom of this file. Each has the form
%     [Nu, Cf] = corr_xxx(Re, c)
% with c a struct of the conditions: c.side ('air' or 'coolant'), c.Pr, c.Tb, c.Tw,
% c.d_h, c.L. Cf must be the FANNING factor (f_Darcy/4), Nu the Nusselt number.
% To add a correlation: copy one of the functions, edit it, and add a row to the
% 'correlations' list below. It then appears in every figure, for air and coolant.
%
% The DNS table starts at the lowest laminar Re it holds; cached results are drawn up
% to its last point (no extrapolated tail). Nu and Cf of the DNS are each one curve.

clear; clc;

%% ---- Correlations to compare (add rows here) ---- %%
% {legend name, function handle}
correlations = {
    'Corr. 1',  @corr_model;
    'Corr. 2 (test 1)',           @corr_test;
    'Corr. 3 (laminar section)',  @corr_lam; 
};

%% ---- Fixed conditions: AIR ---- %%

Pr = 0.715872605452009; %0.703646032163650; %0.715878659971220;
Tb = 2.648278841355512e+02; %3.278969935362533e+02; %2.647969967113464e+02;    % K — bulk air temperature
Tw = 358;                      % K — wall temperature

d_h_air  = 0.002176067743425;            % air hydraulic diameter [m] (same as the model)
L_list   = [0.004 0.05 0.177 inf]; % [d_h_air 0.150378879929866/50 0.150378879929866 1];        % lengths [m] used in place of dx*k_idx
show_dns = true;                        % toggle: draw the DNS curves (air Nu and Cf); false skips the DNS queries too
                                         % (the Nu length factor (1+(d_h/L)^(2/3)) is applied inside thermoturb_cached)

%% ---- Fixed conditions: COOLANT (50/50 EG; no DNS) ---- %%

show_coolant = true;                     % also draw the coolant figures
T_cool    = Tw;                         % K — coolant bulk temperature (property evaluation)
Pr_cool   = mu_EG(T_cool)*cp_EG_50_50(T_cool)/k_EG(T_cool);
r_cool    = 4 - 2*0.1*2/sqrt(3);         % [mm], as r in compare_models
d_h_cool  = d_h_air; %1e-3*r_cool/sqrt(3);         % hydraulic diameter of the equilateral-triangle passage [m]
b_hx      = 1.270442695169008;                         % HX width = coolant flow length [m]   (REPLACE with your case)
L_list_cool = [0.004 0.05 0.5 b_hx inf];                     % lengths [m] for the coolant; e.g. [inf b_hx] to add the real length

show_correlation = true;                 % toggle: draw the correlation curves (air and coolant)

%% ---- Infer Re_tur_actual_min directly from table ---- %%

table_file = fullfile(pwd, 'thermoturb_table.mat');
loaded     = load(table_file, 'tt');
tt         = loaded.tt;

sub_tur           = tt(tt.regime == 1, :);
Re_tur_actual_min = min(sub_tur.Re);
Re_tur_actual_min = max(Re_tur_actual_min, 3000);   % cap at Re_tur_min

fprintf('Turbulent DNS data starts at Re = %.0f\n', Re_tur_actual_min);
fprintf('Table loaded — %d total entries (%d lam, %d tur).\n\n', ...
        height(tt), sum(tt.regime==0), sum(tt.regime==1));

%% ---- Lowest laminar Re in the table (queries typically start around 11-18) ---- %%

Re_lam_min = min(tt.Re(tt.regime == 0));
Re_data_max = max(tt.Re);   % last Re point in the table: the cached results are not drawn beyond it
fprintf('Table Re coverage: %.1f to %.0f\n\n', Re_lam_min, Re_data_max);

%% ---- Re sweep ---- %%

Re_lam_max  = 2300;
Re_tur_min  = 3000;
Re_max_plot = 10000;

Re_arr = unique([...
    linspace(Re_lam_min,        Re_lam_max,        60), ...
    linspace(Re_lam_max,        Re_tur_min,        30), ...
    linspace(Re_tur_min,        Re_tur_actual_min, 30), ...
    linspace(Re_tur_actual_min, Re_max_plot,       50)  ...
]);

n      = numel(Re_arr);
Nu_arr = zeros(n, 1);
Cf_arr = zeros(n, 1);

if show_dns
    fprintf('Querying thermoturb_cached for %d Re values...\n', n);
    fprintf('Pr=%.3f  Tb=%.0fK  Tw=%.0fK\n\n', Pr, Tb, Tw);

    tic
    for i = 1:n
        % L = Inf makes the length factor (1+(d_h/L)^(2/3)) equal 1: the plain DNS result
        [~, Nu_arr(i), Cf_arr(i)] = thermoturb_cached(Re_arr(i), Pr, Tb, Tw, d_h_air, Inf);
    end
    fprintf('Done in %.2fs.\n\n', toc);
else
    Nu_arr(:) = NaN;  Cf_arr(:) = NaN;
end

% DNS Nu with the length factor, as applied inside thermoturb_cached, per L
Nu_dns_L = nan(n, numel(L_list));
if show_dns
    for j = 1:numel(L_list)
        for i = 1:n
            [~, Nu_dns_L(i,j)] = thermoturb_cached(Re_arr(i), Pr, Tb, Tw, d_h_air, L_list(j));
        end
    end
end

% Cached results only up to the table's last point (no extrapolated tail)
beyond = Re_arr(:) > Re_data_max;
Nu_arr(beyond) = NaN;  Cf_arr(beyond) = NaN;  Nu_dns_L(beyond, :) = NaN;

%% ---- Regime masks (DNS Cf) ---- %%

is_lam = Re_arr <= Re_lam_max;
is_tra = Re_arr >  Re_lam_max & Re_arr <  Re_tur_min;
is_ext = Re_arr >= Re_tur_min & Re_arr <  Re_tur_actual_min;
is_tur = Re_arr >= Re_tur_actual_min;

% Handle case where turbulent data extends to or below Re_tur_min
if Re_tur_actual_min <= Re_tur_min
    is_ext = false(n, 1);
    is_tur = Re_arr >= Re_tur_actual_min;
    fprintf('Note: turbulent DNS data extends to Re=%.0f — no extrapolation region.\n\n', ...
            Re_tur_actual_min);
end

%% ---- Correlations: air ---- %%

n_c = size(correlations, 1);
c_air = struct('side', 'air', 'Pr', Pr, 'Tb', Tb, 'Tw', Tw, 'd_h', d_h_air, 'L', NaN);
Nu_corr = nan(n, numel(L_list), n_c);
Cf_corr = nan(n, n_c);
if show_correlation
    for m = 1:n_c
        fn = correlations{m,2};
        for j = 1:numel(L_list)
            c_air.L = L_list(j);
            for i = 1:n
                [Nu_corr(i,j,m), cf_i] = fn(Re_arr(i), c_air);
                if j == 1, Cf_corr(i,m) = cf_i; end     % Cf has no length term
            end
        end
    end
end

%% ---- Correlations: coolant ---- %%

Re_cool = unique([linspace(Re_lam_min, Re_lam_max, 60), linspace(Re_lam_max, Re_tur_min, 30), ...
                  linspace(Re_tur_min, Re_max_plot, 80)]);
n_k = numel(Re_cool);
c_cool = struct('side', 'coolant', 'Pr', Pr_cool, 'Tb', T_cool, 'Tw', T_cool, 'd_h', d_h_cool, 'L', NaN);
Nu_cool = nan(n_k, numel(L_list_cool), n_c);
Cf_cool = nan(n_k, n_c);
if show_coolant && show_correlation
    for m = 1:n_c
        fn = correlations{m,2};
        for j = 1:numel(L_list_cool)
            c_cool.L = L_list_cool(j);
            for i = 1:n_k
                [Nu_cool(i,j,m), cf_i] = fn(Re_cool(i), c_cool);
                if j == 1, Cf_cool(i,m) = cf_i; end
            end
        end
    end
end

%% ---- Colours ---- %%
% One sequential colormap per series: DNS 'Reds', correlation 1 'Blues', 2 'Greens', then
% 'Greys', 'Purples', 'Oranges' (repeating). Within a series, each length gets a shade of
% its map and the LARGEST length is the darkest. Cf (no length term) uses the darkest shade.
cmaps = {'Blues', 'Oranges', 'Greens', 'Purples', 'Greys'};
col_corr   = zeros(n_c*numel(L_list), 3);
col_corr_k = zeros(n_c*numel(L_list_cool), 3);
for m = 1:n_c
    cm = cmaps{mod(m-1, numel(cmaps)) + 1};
    col_corr((m-1)*numel(L_list) + (1:numel(L_list)), :)           = seq_shades(cm, L_list);
    col_corr_k((m-1)*numel(L_list_cool) + (1:numel(L_list_cool)), :) = seq_shades(cm, L_list_cool);
end
col_dns = seq_shades('Reds', L_list);
[~, jmax]   = max(L_list);       % index of the largest length (darkest shade)
[~, jmax_k] = max(L_list_cool);

% Everything the drawing helpers need
D = struct('Re', Re_arr, 'Nu_corr', Nu_corr, 'Cf_corr', Cf_corr, 'names', {correlations(:,1)}, ...
    'L', L_list, 'col_corr', col_corr, 'col_dns', col_dns, 'Nu_dns_L', Nu_dns_L, ...
    'show_dns', show_dns, 'show_corr', show_correlation, 'Cf_dns', Cf_arr, 'jmax', jmax, ...
    'is_lam', is_lam, 'is_tra', is_tra, 'is_ext', is_ext, 'is_tur', is_tur, ...
    'Re_lam_max', Re_lam_max, 'Re_tur_min', Re_tur_min, 'Re_tur_actual_min', Re_tur_actual_min, ...
    'has_dns', true, 'ylim_nu', [0 40], 'ylim_cf', [0 1.4]);
K = struct('Re', Re_cool, 'Nu_corr', Nu_cool, 'Cf_corr', Cf_cool, 'names', {correlations(:,1)}, ...
    'L', L_list_cool, 'col_corr', col_corr_k, 'show_corr', show_correlation, 'jmax', jmax_k, ...
    'Re_lam_max', Re_lam_max, 'Re_tur_min', Re_tur_min, 'has_dns', false, 'ylim_nu', [], 'ylim_cf', [0 1.4]);

air_title  = sprintf('Pr=%.3f  Tb=%.0fK  Tw=%.0fK', Pr, Tb, Tw);
cool_title = sprintf('coolant  Pr=%.2f  T=%.0fK  d_h=%.2f mm', Pr_cool, T_cool, d_h_cool*1e3);

%% ---- AIR: Figure 1 Nu, Figure 2 Cf, Figure 3 combined ---- %%

figure('Name', 'Nu vs Re — air', 'NumberTitle', 'off');
ax = axes; hold(ax, 'on');
draw_nu(ax, D, false);
title(ax, ['Nu vs Re  |  ' air_title]);  legend(ax, 'Location', 'northwest');

figure('Name', 'Cf vs Re — air', 'NumberTitle', 'off');
ax = axes; hold(ax, 'on');
draw_cf(ax, D, false);
title(ax, ['Cf vs Re  |  ' air_title]);  legend(ax, 'Location', 'northeast');

figure('Name', 'Nu and Cf vs Re — air', 'NumberTitle', 'off');
ax1 = subplot(2,1,1); hold(ax1, 'on');
draw_nu(ax1, D, true);
title(ax1, ['Nu and Cf vs Re  |  ' air_title]);  legend(ax1, 'Location', 'northwest', 'FontSize', 8);
ax2 = subplot(2,1,2); hold(ax2, 'on');
draw_cf(ax2, D, true);
legend(ax2, 'Location', 'northeast', 'FontSize', 8);

%% ---- COOLANT: Figure 4 Nu, Figure 5 Cf, Figure 6 combined ---- %%

if show_coolant
    figure('Name', 'Nu vs Re — coolant', 'NumberTitle', 'off');
    ax = axes; hold(ax, 'on');
    draw_nu(ax, K, false);
    title(ax, ['Nu vs Re  |  ' cool_title]);  legend(ax, 'Location', 'northwest');

    figure('Name', 'Cf vs Re — coolant', 'NumberTitle', 'off');
    ax = axes; hold(ax, 'on');
    draw_cf(ax, K, false);
    title(ax, ['Cf vs Re  |  ' cool_title]);  legend(ax, 'Location', 'northeast');

    figure('Name', 'Nu and Cf vs Re — coolant', 'NumberTitle', 'off');
    ax1 = subplot(2,1,1); hold(ax1, 'on');
    draw_nu(ax1, K, true);
    title(ax1, ['Nu and Cf vs Re  |  ' cool_title]);  legend(ax1, 'Location', 'northwest', 'FontSize', 8);
    ax2 = subplot(2,1,2); hold(ax2, 'on');
    draw_cf(ax2, K, true);
    legend(ax2, 'Location', 'northeast', 'FontSize', 8);
end

%% ---- Print summary table (air DNS) ---- %%

if show_dns

fprintf('%-10s %-12s %-12s %-12s\n', 'Re', 'Nu (DNS)', 'Cf (DNS)', 'Regime');
fprintf('%s\n', repmat('-', 1, 50));
step = max(1, floor(n/15));
for i = 1:step:n
    if is_lam(i),      reg = 'laminar';
    elseif is_tra(i),  reg = 'transition';
    elseif any(is_ext) && is_ext(i), reg = 'tur-extrap';
    else,              reg = 'turbulent';
    end
    fprintf('%-10.0f %-12.4f %-12.4e %-12s\n', Re_arr(i), Nu_arr(i), Cf_arr(i), reg);
end
end

%% ======================= CORRELATIONS (edit / add here) ======================= %%
% Each function: [Nu, Cf] = corr_xxx(Re, c); Cf = Fanning friction factor.
% c.side = 'air' or 'coolant'; c.Pr, c.Tb, c.Tw [K], c.d_h [m], c.L [m] (flow length; Inf = fully developed).

%% ---- Correlation 1: as in HX_design1_disc ---- %%
% Air:     Nu = 3.66 (Re < 2300), linear blend 2300-3000, Gnielinski above, with
%          f = Konakov (Fanning, f/2 in Gnielinski), entrance factor (1 + 1/3 (d_h/L)^(2/3)),
%          temperature factor K_ratio^0.45 with K_ratio = Tb/Tw. Cf = f.
% Coolant: Nu = 4.36 (Re < 2300), Nu = 0.023 Re^0.8 Pr^0.4 above; f_Darcy from
%          (-1.8 log10((0.0015/d_h)^1.11 + 6.9/Re))^-2, Cf = f_Darcy/4.
function [Nu, Cf] = corr_model(Re, c)
sr = 2e-06;  
if strcmp(c.side, 'air')
    f = 0.25*((1.8*log10(Re)) - 1.5)^(-2);          % Konakov, Fanning
    K_ratio = (c.Tb/c.Tw)^0.45;
    Nu_turb = (f/2)*(Re-1000)*c.Pr*K_ratio*(1 + 1/3*(c.d_h/c.L)^(2/3)) / (1 + 12.7*sqrt(f/2)*(c.Pr^(2/3) - 1));
    if Re < 2300
        Nu = 3.66;
    elseif Re < 3000
        w  = (Re - 2300)/(3000 - 2300);
        Nu = (1-w)*3.66 + w*Nu_turb;
    else
        Nu = Nu_turb;
    end
    Cf = f;
else
    f_D = (1/(-1.8*log10((sr/c.d_h)^1.11 + (6.9/Re))))^2;
    if Re < 2300
        Nu = 4.36;
    else
        Nu = 0.023*Re^0.8*c.Pr^0.4;
    end
    Cf = f_D/4;
end
end

%% ---- Correlation 2: test set (lumped HX_design1_v2 style) ---- %%
% f_Darcy: Konakov below Re 1800, Gerl above 3000, cosine blend between. Cf = f_Darcy/4.
% Air:     Nu laminar = 4.364 + 0.086 (Re Pr d_h/L)^1.33 / (1 + 0.1 Pr (Re d_h/L)^0.83);
%          Gnielinski with f/8, sqrt(f/8), K^0.45, K = Tb/Tw, blended linearly 2300-3000.
% Coolant: same laminar Nu with L; Gnielinski with f/8 and sqrt(f/2) (as in the file), K = 1.
function [Nu, Cf] = corr_test(Re, c)
sr = 2e-06;                                          % m, aluminium surface roughness
Re_low = 1800;  Re_high = 3000;
f_lam  = ((1.8*log10(Re)) - 1.5)^(-2);               % Konakov (Darcy)
f_turb = (1/(-2*log10(2.7*log10(Re)^1.2/Re + (sr/c.d_h)/3.71)))^2;   % Gerl (Darcy)
if Re <= Re_low
    f = f_lam;
elseif Re >= Re_high
    f = f_turb;
else
    x = (Re - Re_low)/(Re_high - Re_low);
    w = 0.5*(1 - cos(pi*x));
    f = (1-w)*f_lam + w*f_turb;
end
Cf = f/4;
if strcmp(c.side, 'air')
    K = (c.Tb/c.Tw)^0.45;                                   % T_mean_c / T_mean_h
else
    K = 1;
end
Nu_lam = 2.49 + (0.0677*(Re*c.Pr*c.d_h/c.L)^1.33)/(1 + 0.1*c.Pr*(c.d_h*Re/c.L)^0.83);
Nu_turb = (f/8)*(Re-1000)*c.Pr*(1 + 1/3*(c.d_h/c.L)^(2/3))*K / (1 + 12.7*sqrt(f/8)*(c.Pr^(2/3) - 1));
if Re < 2300
    Nu = Nu_lam;
elseif Re < 3000
    w  = (Re - 2300)/(3000 - 2300);
    Nu = (1-w)*Nu_lam + w*Nu_turb;
else
    Nu = Nu_turb;
end
end

%% ---- Correlation 3: test set 2 (inclusion of laminar entry) ---- %%
% f_Darcy: Konakov below Re 1800, Gerl above 3000, cosine blend between. Cf = f_Darcy/4.
% Air:     Nu laminar = 4.364 + 0.086 (Re Pr d_h/L)^1.33 / (1 + 0.1 Pr (Re d_h/L)^0.83);
%          Gnielinski with f/8, sqrt(f/8), K^0.45, K = Tb/Tw, blended linearly 2300-3000.
% Coolant: same laminar Nu with L; Gnielinski with f/8 and sqrt(f/2) (as in the file), K = 1.
function [Nu, Cf] = corr_lam(Re, c)
sr = 2e-06;                                          % m, aluminium surface roughness
Re_low = 1800;  Re_high = 3000;

Nu_x_1 = 2.49;
Nu_x_2 = 1.077*c.Pr^(1/3)*sqrt(Re*c.d_h/c.L);
Nu_x_3 = 0.5*(2/(1+22*c.Pr))^(1/6) * sqrt(Re*c.Pr*c.d_h/c.L);

% Nu_x_lam = (Nu_x_1^3 + 0.7^3 + (Nu_x_2-0.7)^3 + Nu_x_3^3)^(1/3);

f_lam  = ((1.8*log10(Re)) - 1.5)^(-2);               % Konakov (Darcy)
f_turb = (1/(-2*log10(2.7*log10(Re)^1.2/Re + (sr/c.d_h)/3.71)))^2;   % Gerl (Darcy)

Ref_lam_fd = 40/3;
x_star = c.L/Re/c.Pr/c.d_h;
f_lam = (Ref_lam_fd + 25/sqrt(x_star))/Re;
Nu_x_lam = 0.4915*(x_star)^-0.3286;

f_lam = 53.3/Re;
if Re <= 1800
    f = f_lam;
elseif Re >= 3000
    f = f_turb;
    % f = 0.0425 / Re^0.2;
else
    wb = 0.5*(1-cos(pi*(Re-1800)/(3000-1800)));   % cosine smoothing
    f  = (1-wb)*f_lam + wb*f_turb;
end
Cf  = f/4;
% if Re <= Re_low
%     f = f_lam;
% elseif Re >= Re_high
%     f = f_turb;
% else
%     x = (Re - Re_low)/(Re_high - Re_low);
%     w = 0.5*(1 - cos(pi*x));
%     f = (1-w)*f_lam + w*f_turb;
% end
% Cf = f/4;
if strcmp(c.side, 'air')
    K = (c.Tb/c.Tw)^0.45;                                   % T_mean_c / T_mean_h
else
    K = 1;
end
% Nu_lam = 4.364 + (0.086*(Re*c.Pr*c.d_h/c.L)^1.33)/(1 + 0.1*c.Pr*(c.d_h*Re/c.L)^0.83);
Nu_lam = Nu_x_lam;
Nu_turb = (f/8)*(Re-1000)*c.Pr*(1 + 1/3*(c.d_h/c.L)^(2/3))*K / (1 + 12.7*sqrt(f/8)*(c.Pr^(2/3) - 1));
if Re < 2300
    Nu = Nu_lam;
elseif Re < 3000
    w  = (Re - 2300)/(3000 - 2300);
    Nu = (1-w)*Nu_lam + w*Nu_turb;
else
    Nu = Nu_turb;
end
end

%% ======================= PLOTTING HELPERS ======================= %%

%% ---- Local function: line style by length ---- %%
function ls = ls_of(L)
% Fully developed (L = Inf): solid. Any finite length: striped (dash-dot).
if isinf(L), ls = '-'; else, ls = '-.'; end
end

%% ---- Local function: vertical regime lines ---- %%
function regime_lines(S, short)
xline(S.Re_lam_max, '--k', 'Re=2300', 'LabelHorizontalAlignment', 'left',  'LineWidth', 1, 'HandleVisibility', 'off')
xline(S.Re_tur_min, '--k', 'Re=3000', 'LabelHorizontalAlignment', 'right', 'LineWidth', 1, 'HandleVisibility', 'off')
if S.has_dns && S.show_dns && S.Re_tur_actual_min > S.Re_tur_min
    if short, lab = sprintf('Re=%.0f', S.Re_tur_actual_min); else, lab = sprintf('Re=%.0f (data)', S.Re_tur_actual_min); end
    xline(S.Re_tur_actual_min, '--k', lab, 'LabelHorizontalAlignment', 'right', 'LineWidth', 1, 'HandleVisibility', 'off')
end
end

%% ---- Local function: Nu vs Re panel ---- %%
function draw_nu(ax, S, short)
regime_lines(S, short);
n_L = numel(S.L);
if S.show_corr
    for m = 1:numel(S.names)
        for j = 1:n_L
            plot(ax, S.Re, S.Nu_corr(:,j,m), ls_of(S.L(j)), 'Color', S.col_corr((m-1)*n_L + j, :), 'LineWidth', 2, ...
                 'DisplayName', sprintf('%s, L = %g m', S.names{m}, S.L(j)))
        end
    end
end
if S.has_dns && S.show_dns
    for j = 1:n_L
        if short, lab = 'DNS + factor'; else, lab = 'DNS with length factor'; end
        plot(ax, S.Re, S.Nu_dns_L(:,j), ls_of(S.L(j)), 'Color', S.col_dns(j,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('%s, L = %g m', lab, S.L(j)))
    end
end
xlabel(ax, 'Re_{bulk}')
ylabel(ax, 'Nusselt number [-]')
grid(ax, 'on')
if ~isempty(S.ylim_nu), ylim(ax, S.ylim_nu), end
end

%% ---- Local function: Cf vs Re panel ---- %%
function draw_cf(ax, S, short)
if S.has_dns && S.show_dns
    % One DNS curve over the whole Re range (cached results up to the table's last point)
    plot(ax, S.Re, S.Cf_dns, '-', 'Color', S.col_dns(S.jmax,:), 'LineWidth', 2, 'DisplayName', 'DNS (thermoturb\_cached)')
end
regime_lines(S, short);
if S.show_corr
    n_L = numel(S.L);
    for m = 1:numel(S.names)
        plot(ax, S.Re, S.Cf_corr(:,m), '-', 'Color', S.col_corr((m-1)*n_L + S.jmax, :), 'LineWidth', 2, ...
             'DisplayName', S.names{m})
    end
end
xlabel(ax, 'Re_{bulk}')
ylabel(ax, 'C_f (Fanning) [-]')
grid(ax, 'on')
if ~isempty(S.ylim_cf), ylim(ax, S.ylim_cf), end
end

%% ---- Local function: shades of a sequential colormap, one per length ---- %%
% The largest length gets the darkest shade, the smallest a lighter one (single length: dark).
function C = seq_shades(name, L)
switch name
    case 'Reds'
        A = [255 245 240; 254 224 210; 252 187 161; 252 146 114; 251 106 74; 239 59 44; 203 24 29; 165 15 21; 103 0 13];
    case 'Blues'
        A = [247 251 255; 222 235 247; 198 219 239; 158 202 225; 107 174 214; 66 146 198; 33 113 181; 8 81 156; 8 48 107];
    case 'Greys'
        A = [255 255 255; 240 240 240; 217 217 217; 189 189 189; 150 150 150; 115 115 115; 82 82 82; 37 37 37; 0 0 0];
    case 'Greens'
        A = [247 252 245; 229 245 224; 199 233 192; 161 217 155; 116 196 118; 65 171 93; 35 139 69; 0 109 44; 0 68 27];
    case 'Purples'
        A = [252 251 253; 239 237 245; 218 218 235; 188 189 220; 158 154 200; 128 125 186; 106 81 163; 84 39 143; 63 0 125];
    case 'Oranges'
        A = [255 245 235; 254 230 206; 253 208 162; 253 174 107; 253 141 60; 241 105 19; 217 72 1; 166 54 3; 127 39 4];
    otherwise
        error('seq_shades:unknown', 'Unknown colormap "%s".', name);
end
nL = numel(L);
[~, ord] = sort(L(:));
rk = zeros(nL, 1);  rk(ord) = 1:nL;                % rank 1 = smallest length, nL = largest
if nL == 1
    t = 0.8;
else
    t = 0.4 + 0.55*(rk - 1)/(nL - 1);               % 0.4 (lightest used) ... 0.95 (darkest)
end
x = linspace(0, 1, size(A, 1));
C = interp1(x, A/255, t);
end