%% test_thermoturb_sweep.m
% Tests thermoturb_query_sweep for one (Pr, Tb, Tw) combination
% and plots the resulting Nu and Cf vs Re curves.

clear; clc;

%% ---- Test inputs ---- %%

Pr = 0.703646032163650;
Tb = 3.278969935362533e+02;    % K — bulk air temperature
Tw = 358;                      % K — wall temperature

%% ---- Correlation comparison (from calculate_air_side in HX_design1_disc) ---- %%

show_correlation = true;                 % toggle: also draw the correlation Nu curves
d_h_air  = 0.002209401076759;            % air hydraulic diameter [m] (same as the model)
L_list   = [inf];                        % lengths [m] for the entrance factor (Inf = none)
K_ratio  = (Tb / Tw)^0.45;                      % model uses T_hot_ref / T_in, with T_hot_ref = Tw here
col_corr = autumn(numel(L_list) + 2);    % correlation: red-yellow, one colour per length (skip the extremes)
col_corr = col_corr(2:end-1, :);

flow_model = 'laminar';
Re_min     = 10;
Re_max     = 5000;

fprintf('Testing thermoturb_query_sweep:\n');
fprintf('  Pr = %.3f\n', Pr);
fprintf('  Tb = %.0f K\n', Tb);
fprintf('  Tw = %.0f K\n', Tw);
fprintf('  Flow model: %s\n', flow_model);
fprintf('  Re range: [%d - %d]\n\n', Re_min, Re_max);

%% ---- Run laminar sweep ---- %%

tic
sweep = thermoturb_query_sweep(Re_min, Re_max, Pr, Tb, Tw, flow_model);
elapsed = toc;

fprintf('\nQuery completed in %.2f s\n', elapsed);
fprintf('Returned %d Re points\n', height(sweep));

if height(sweep) > 10
    fprintf('Sweep mode confirmed.\n\n');
else
    fprintf('WARNING: few points returned — check generateCoefficients field.\n\n');
end

%% ---- Print summary ---- %%

fprintf('%-12s %-12s %-12s %-12s %-12s %-12s\n', ...
        'Re_bulk', 'Nu', 'Cf', 'Stanton', 'Retau', 'Retau_cp');
fprintf('%s\n', repmat('-', 1, 72));
step = max(1, floor(height(sweep)/10));   % print ~10 evenly spaced rows
for i = 1:step:height(sweep)
    fprintf('%-12.1f %-12.4f %-12.6f %-12.6f %-12.4f %-12.4f\n', ...
            sweep.Re_bulk(i), sweep.Nu(i), sweep.Cf(i), ...
            sweep.Stanton(i), sweep.Retau(i), sweep.Retau_cp(i));
end

%% ---- Plots ---- %%

figure('Name', 'thermoturb sweep — Nu and Cf vs Re', 'NumberTitle', 'off');

subplot(2,1,1)
plot(sweep.Re_bulk, sweep.Nu, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'DNS (laminar sweep)')
hold on
if show_correlation
    for j = 1:numel(L_list)
        plot(sweep.Re_bulk, nu_correlation(sweep.Re_bulk, Pr, K_ratio, d_h_air, L_list(j)), ls_of(L_list(j)), ...
             'Color', col_corr(j,:), 'LineWidth', 1.3, 'DisplayName', sprintf('Correlation, L = %g m', L_list(j)))
    end
    legend('Location', 'best')
end
xlabel('Re_{bulk}')
ylabel('Nusselt number [-]')
title(sprintf('Nu vs Re  |  %s  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', ...
              flow_model, Pr, Tb, Tw))
grid on
ylim([0 40])

subplot(2,1,2)
plot(sweep.Re_bulk, sweep.Cf, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'DNS (laminar sweep)')
hold on
if show_correlation
    plot(sweep.Re_bulk, cf_correlation(sweep.Re_bulk), '-', 'Color', col_corr(1,:), 'LineWidth', 1.3, ...
         'DisplayName', 'Correlation (Konakov, as in the model)')
    legend('Location', 'best')
end
xlabel('Re_{bulk}')
ylabel('C_f (Fanning) [-]')
title(sprintf('Cf vs Re  |  %s  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', ...
              flow_model, Pr, Tb, Tw))
grid on
ylim([0 1.4])

%% ---- Also test turbulent regime ---- %%

fprintf('\nNow testing turbulent regime Re=[5400-10000]...\n');

tic
sweep_tur = thermoturb_query_sweep(5400, 10000, Pr, Tb, Tw, 'turbulent');
elapsed_tur = toc;

fprintf('Query completed in %.2f s\n', elapsed_tur);
fprintf('Returned %d Re points\n\n', height(sweep_tur));

%% ---- Print summary ---- %%

fprintf('%-12s %-12s %-12s %-12s %-12s %-12s\n', ...
    'Re_bulk', 'Nu', 'Cf', 'Stanton', 'Retau', 'Retau_cp');
fprintf('%s\n', repmat('-', 1, 72));
step = max(1, floor(height(sweep_tur)/10));   % print ~10 evenly spaced rows
for i = 1:step:height(sweep_tur)
    fprintf('%-12.1f %-12.4f %-12.6f %-12.6f %-12.4f %-12.4f\n', ...
        sweep_tur.Re_bulk(i), sweep_tur.Nu(i), sweep_tur.Cf(i), ...
        sweep_tur.Stanton(i), sweep_tur.Retau(i), sweep_tur.Retau_cp(i));
end

%% ---- Plots ---- %%

figure('Name', 'thermoturb sweep — Nu and Cf vs Re', 'NumberTitle', 'off');

subplot(2,1,1)
plot(sweep_tur.Re_bulk, sweep_tur.Nu, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'DNS (turbulent sweep)')
hold on
if show_correlation
    for j = 1:numel(L_list)
        plot(sweep_tur.Re_bulk, nu_correlation(sweep_tur.Re_bulk, Pr, K_ratio, d_h_air, L_list(j)), ls_of(L_list(j)), ...
             'Color', col_corr(j,:), 'LineWidth', 1.3, 'DisplayName', sprintf('Correlation, L = %g m', L_list(j)))
    end
    legend('Location', 'best')
end
xlabel('Re_{bulk}')
ylabel('Nusselt number [-]')
title(sprintf('Nu vs Re  |  turbulent  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', ...
    Pr, Tb, Tw))
grid on
ylim([0 40])

subplot(2,1,2)
plot(sweep_tur.Re_bulk, sweep_tur.Cf, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'DNS (turbulent sweep)')
hold on
if show_correlation
    plot(sweep_tur.Re_bulk, cf_correlation(sweep_tur.Re_bulk), '-', 'Color', col_corr(1,:), 'LineWidth', 1.3, ...
         'DisplayName', 'Correlation (Konakov, as in the model)')
    legend('Location', 'best')
end
xlabel('Re_{bulk}')
ylabel('C_f (Fanning) [-]')
title(sprintf('Cf vs Re  |  turbulent  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', ...
    Pr, Tb, Tw))
grid on
ylim([0 1.4])



%% ---- Plot combined regimes --- %%

figure('Name', 'thermoturb sweep — both regimes', 'NumberTitle', 'off');

subplot(2,1,1)
hold on
plot(sweep.Re_bulk,     sweep.Nu,     'b-o', 'LineWidth', 1.5, ...
     'MarkerSize', 4, 'DisplayName', 'Laminar')
plot(sweep_tur.Re_bulk, sweep_tur.Nu, 'r-o', 'LineWidth', 1.5, ...
     'MarkerSize', 4, 'DisplayName', 'Turbulent')
if show_correlation
    Re_corr = [sweep.Re_bulk(:); sweep_tur.Re_bulk(:)];
    Re_corr = sort(Re_corr);
    for j = 1:numel(L_list)
        plot(Re_corr, nu_correlation(Re_corr, Pr, K_ratio, d_h_air, L_list(j)), ls_of(L_list(j)), ...
             'Color', col_corr(j,:), 'LineWidth', 1.3, 'DisplayName', sprintf('Correlation, L = %g m', L_list(j)))
    end
end
xlabel('Re_{bulk}')
ylabel('Nusselt number [-]')
title(sprintf('Nu vs Re  |  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', Pr, Tb, Tw))
legend('Location', 'best')
grid on
ylim([0 40])

subplot(2,1,2)
hold on
plot(sweep.Re_bulk,     sweep.Cf,     'b-o', 'LineWidth', 1.5, ...
     'MarkerSize', 4, 'DisplayName', 'Laminar')
plot(sweep_tur.Re_bulk, sweep_tur.Cf, 'r-o', 'LineWidth', 1.5, ...
     'MarkerSize', 4, 'DisplayName', 'Turbulent')
if show_correlation
    plot(Re_corr, cf_correlation(Re_corr), '-', 'Color', col_corr(1,:), 'LineWidth', 1.3, ...
         'DisplayName', 'Correlation (Konakov, as in the model)')
end
xlabel('Re_{bulk}')
ylabel('C_f (Fanning) [-]')
title(sprintf('Cf vs Re  |  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', Pr, Tb, Tw))
legend('Location', 'best')
grid on
ylim([0 1.4])

%% ---- Local function: correlation Nu(Re), as in calculate_air_side (no DNS) ---- %%

function Nu = nu_correlation(Re, Pr, K_ratio, d_h, L)
% Laminar 3.66 below Re 2300, Konakov/Gnielinski-type turbulent form above 3000,
% linear blend in between. L = Inf removes the entrance (length) factor.
Nu = zeros(size(Re));
for i = 1:numel(Re)
    Re_in = Re(i);
    f_in  = 0.25*((1.8*log10(Re_in))-1.5)^(-2);   % Konakov et al, as in the model
    Nu_turb = (f_in/2)*(Re_in-1000)*Pr*K_ratio* ...
              (1 + 1/3*(d_h/L)^(2/3)) / ...
              (1+12.7*sqrt(f_in/2)*(Pr^(2/3)-1));
    if Re_in < 2300
        Nu(i) = 3.66;
    elseif Re_in < 3000
        w = (Re_in-2300)/(3000-2300);
        Nu(i) = (1-w)*3.66 + w*Nu_turb;
    else
        Nu(i) = Nu_turb;
    end
end
end

function ls = ls_of(L)
% Fully developed (L = Inf): solid. Any finite length: striped (dash-dot).
if isinf(L), ls = '-'; else, ls = '-.'; end
end

function Cf = cf_correlation(Re)
% Friction factor the model uses at every Re (Konakov et al, Fanning); no length term.
Cf = 0.25*((1.8*log10(Re))-1.5).^(-2);
end