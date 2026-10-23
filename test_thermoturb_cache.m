%% test_thermoturb_cache.m
% Plots Nu and Cf from thermoturb_cached over the full Re range,
% showing laminar, transition, extrapolated turbulent, and DNS turbulent
% regimes continuously. Infers turbulent data boundary from table directly.

clear; clc;

%% ---- Fixed conditions ---- %%

Pr = 0.715872605452009; %0.703646032163650; %0.715878659971220;
Tb = 2.648278841355512e+02; %3.278969935362533e+02; %2.647969967113464e+02;    % K — bulk air temperature
Tw = 358;                      % K — wall temperature

%% ---- Correlation comparison (from calculate_air_side in HX_design1_disc) ---- %%

show_correlation = true;                 % toggle: also draw the correlation Nu curves
d_h_air  = 0.002209401076759;            % air hydraulic diameter [m] (same as the model)
% f_in is computed from Re at each point (Konakov, as in the model), see below
L_list   = [inf]; % [d_h_air 0.150378879929866/50 0.150378879929866 1];        % lengths [m] used in place of dx*k_idx
K_ratio  = (Tb / Tw)^(0.45);                      % model uses T_hot_ref / T_in, with T_hot_ref = Tw here
show_dns_scaled = true;                  % keep true: the Nu plots show only this DNS curve (not the plain cached one). Toggle: also draw the DNS Nu from thermoturb_cached(..., d_h, L) for each L in L_list
                                         % (the length factor (1+(d_h/L)^(2/3)) is applied inside thermoturb_cached; no K_ratio)

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

fprintf('Querying thermoturb_cached for %d Re values...\n', n);
fprintf('Pr=%.3f  Tb=%.0fK  Tw=%.0fK\n\n', Pr, Tb, Tw);

tic
for i = 1:n
    % L = Inf makes the length factor (1+(d_h/L)^(2/3)) equal 1: the plain DNS result
    [~, Nu_arr(i), Cf_arr(i)] = thermoturb_cached(Re_arr(i), Pr, Tb, Tw, d_h_air, Inf);
end
fprintf('Done in %.2fs.\n\n', toc);

% DNS Nu with the length factor, as applied inside thermoturb_cached, per L
Nu_dns_L = nan(n, numel(L_list));
if show_dns_scaled
    for j = 1:numel(L_list)
        for i = 1:n
            [~, Nu_dns_L(i,j)] = thermoturb_cached(Re_arr(i), Pr, Tb, Tw, d_h_air, L_list(j));
        end
    end
end

% Cached results only up to the table's last point (no extrapolated tail)
beyond = Re_arr(:) > Re_data_max;
Nu_arr(beyond) = NaN;  Cf_arr(beyond) = NaN;  Nu_dns_L(beyond, :) = NaN;

%% ---- Regime masks ---- %%

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

%% ---- Correlation Nu(Re) for each length ---- %%

Nu_corr = zeros(n, numel(L_list));
if show_correlation
    for j = 1:numel(L_list)
        for i = 1:n
            Re_in = Re_arr(i);
            f_in  = 0.25*((1.8*log10(Re_in))-1.5)^(-2);   % Konakov et al, as in the model
            Nu_turb = (f_in/2)*(Re_in)*Pr*K_ratio* ...
                      (1 + 1/3*(d_h_air/L_list(j))^(2/3)) / ... % (1+(d_h_air/L_list(j))^(2/3))*K_ratio / ...
                      (1+12.7*sqrt(f_in/2)*(Pr^(2/3)-1));
            if Re_in < 2300
                Nu_corr(i,j) = 3.66;
            elseif Re_in < 3000
                w = (Re_in-2300)/(3000-2300);
                Nu_corr(i,j) = (1-w)*3.66 + w*Nu_turb;
            else
                Nu_corr(i,j) = Nu_turb;
            end
        end
    end
end
Cf_corr = 0.25*((1.8*log10(Re_arr))-1.5).^(-2);   % friction factor the model uses at every Re (Konakov, Fanning); no length term
col_corr = autumn(numel(L_list) + 2);    % correlation: red-yellow, one colour per length (skip the extremes)
col_corr = col_corr(2:end-1, :);
col_dns  = winter(numel(L_list) + 2);    % DNS with length factor: blue-green, one colour per length
col_dns  = col_dns(2:end-1, :);

%% ---- Colour scheme ---- %%

col_lam = [0.122 0.471 0.706];   % blue   — laminar DNS
col_tra = [0.839 0.153 0.157];   % red    — transition cosine blend
col_ext = [0.580 0.404 0.741];   % purple — turbulent extrapolated
col_tur = [0.173 0.627 0.173];   % green  — turbulent DNS
lw = 2;

%% ---- Figure 1: Nu vs Re ---- %%

figure('Name', 'Nu vs Re — thermoturb_cached', 'NumberTitle', 'off');
hold on

xline(Re_lam_max, '--k', 'Re=2300', 'LabelHorizontalAlignment', 'left',  'LineWidth', 1, 'HandleVisibility', 'off')
xline(Re_tur_min, '--k', 'Re=3000', 'LabelHorizontalAlignment', 'right', 'LineWidth', 1, 'HandleVisibility', 'off')
if Re_tur_actual_min > Re_tur_min
    xline(Re_tur_actual_min, '--k', sprintf('Re=%.0f (data)', Re_tur_actual_min), ...
          'LabelHorizontalAlignment', 'right', 'LineWidth', 1, 'HandleVisibility', 'off')
end

if show_correlation
    for j = 1:numel(L_list)
        plot(Re_arr, Nu_corr(:,j), ls_of(L_list(j)), 'Color', col_corr(j,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('Correlation, L = %g m', L_list(j)))
    end
end

if show_dns_scaled
    for j = 1:numel(L_list)
        plot(Re_arr, Nu_dns_L(:,j), ls_of(L_list(j)), 'Color', col_dns(j,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('DNS with length factor, L = %g m', L_list(j)))
    end
end

xlabel('Re_{bulk}')
ylabel('Nusselt number [-]')
title(sprintf('Nu vs Re  |  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', Pr, Tb, Tw))
legend('Location', 'northwest')
grid on
ylim([0 40])

%% ---- Figure 2: Cf vs Re ---- %%

figure('Name', 'Cf vs Re — thermoturb_cached', 'NumberTitle', 'off');
hold on

plot(Re_arr(is_lam), Cf_arr(is_lam), '-',  'Color', col_lam, 'LineWidth', lw, ...
     'DisplayName', 'Laminar (DNS)')
plot(Re_arr(is_tra), Cf_arr(is_tra), '--', 'Color', col_tra, 'LineWidth', lw, ...
     'DisplayName', 'Transition (cosine blend)')
if any(is_ext)
    plot(Re_arr(is_ext), Cf_arr(is_ext), ':', 'Color', col_ext, 'LineWidth', lw, ...
         'DisplayName', sprintf('Turbulent (extrap. from Re=%.0f)', Re_tur_actual_min))
end
plot(Re_arr(is_tur), Cf_arr(is_tur), '-',  'Color', col_tur, 'LineWidth', lw, ...
     'DisplayName', sprintf('Turbulent DNS (Re\\geq%.0f)', Re_tur_actual_min))

xline(Re_lam_max, '--k', 'Re=2300', 'LabelHorizontalAlignment', 'left',  'LineWidth', 1, 'HandleVisibility', 'off')
xline(Re_tur_min, '--k', 'Re=3000', 'LabelHorizontalAlignment', 'right', 'LineWidth', 1, 'HandleVisibility', 'off')
if Re_tur_actual_min > Re_tur_min
    xline(Re_tur_actual_min, '--k', sprintf('Re=%.0f (data)', Re_tur_actual_min), ...
          'LabelHorizontalAlignment', 'right', 'LineWidth', 1, 'HandleVisibility', 'off')
end

if show_correlation
    plot(Re_arr, Cf_corr, '-', 'Color', col_corr(1,:), 'LineWidth', 2, ...
         'DisplayName', 'Correlation (Konakov, as in the model)')
end

xlabel('Re_{bulk}')
ylabel('C_f (Fanning) [-]')
title(sprintf('Cf vs Re  |  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', Pr, Tb, Tw))
legend('Location', 'northeast')
grid on
ylim([0 1.4])

%% ---- Figure 3: Combined 2x1 panel ---- %%

figure('Name', 'Nu and Cf vs Re — thermoturb_cached', 'NumberTitle', 'off');

subplot(2,1,1)
hold on
xline(Re_lam_max, '--k', 'Re=2300', 'LabelHorizontalAlignment', 'left',  'LineWidth', 1, 'HandleVisibility', 'off')
xline(Re_tur_min, '--k', 'Re=3000', 'LabelHorizontalAlignment', 'right', 'LineWidth', 1, 'HandleVisibility', 'off')
if Re_tur_actual_min > Re_tur_min
    xline(Re_tur_actual_min, '--k', sprintf('Re=%.0f', Re_tur_actual_min), ...
          'LabelHorizontalAlignment', 'right', 'LineWidth', 1, 'HandleVisibility', 'off')
end
if show_correlation
    for j = 1:numel(L_list)
        plot(Re_arr, Nu_corr(:,j), ls_of(L_list(j)), 'Color', col_corr(j,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('Corr., L = %g m', L_list(j)))
    end
end
if show_dns_scaled
    for j = 1:numel(L_list)
        plot(Re_arr, Nu_dns_L(:,j), ls_of(L_list(j)), 'Color', col_dns(j,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('DNS + factor, L = %g m', L_list(j)))
    end
end
ylabel('Nusselt number [-]')
title(sprintf('Nu and Cf vs Re  |  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', Pr, Tb, Tw))
legend('Location', 'northwest', 'FontSize', 8)
grid on
ylim([0 40])

subplot(2,1,2)
hold on
plot(Re_arr(is_lam), Cf_arr(is_lam), '-',  'Color', col_lam, 'LineWidth', lw, ...
     'DisplayName', 'Laminar (DNS)')
plot(Re_arr(is_tra), Cf_arr(is_tra), '--', 'Color', col_tra, 'LineWidth', lw, ...
     'DisplayName', 'Transition')
if any(is_ext)
    plot(Re_arr(is_ext), Cf_arr(is_ext), ':', 'Color', col_ext, 'LineWidth', lw, ...
         'DisplayName', 'Turbulent (extrap.)')
end
plot(Re_arr(is_tur), Cf_arr(is_tur), '-',  'Color', col_tur, 'LineWidth', lw, ...
     'DisplayName', 'Turbulent (DNS)')
xline(Re_lam_max, '--k', 'Re=2300', 'LabelHorizontalAlignment', 'left',  'LineWidth', 1, 'HandleVisibility', 'off')
xline(Re_tur_min, '--k', 'Re=3000', 'LabelHorizontalAlignment', 'right', 'LineWidth', 1, 'HandleVisibility', 'off')
if Re_tur_actual_min > Re_tur_min
    xline(Re_tur_actual_min, '--k', sprintf('Re=%.0f', Re_tur_actual_min), ...
          'LabelHorizontalAlignment', 'right', 'LineWidth', 1, 'HandleVisibility', 'off')
end
if show_correlation
    plot(Re_arr, Cf_corr, '-', 'Color', col_corr(1,:), 'LineWidth', 2, ...
         'DisplayName', 'Corr. (Konakov)')
end
xlabel('Re_{bulk}')
ylabel('C_f (Fanning) [-]')
legend('Location', 'northeast', 'FontSize', 8)
grid on
ylim([0 1.4])

%% ---- Print summary table ---- %%

fprintf('%-10s %-12s %-12s %-12s\n', 'Re', 'Nu', 'Cf', 'Regime');
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

%% ---- Local function: line style by length ---- %%

function ls = ls_of(L)
% Fully developed (L = Inf): solid. Any finite length: striped (dash-dot).
if isinf(L), ls = '-'; else, ls = '-.'; end
end