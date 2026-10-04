%% test_thermoturb_sweep.m
% Tests thermoturb_query_sweep for one (Pr, Tb, Tw) combination
% and plots the resulting Nu and Cf vs Re curves.

clear; clc;

%% ---- Test inputs ---- %%

Pr         = 0.7;
Tb         = 280;    % K — bulk air temperature
Tw         = 350;    % K — wall temperature
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
plot(sweep.Re_bulk, sweep.Nu, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4)
xlabel('Re_{bulk}')
ylabel('Nusselt number [-]')
title(sprintf('Nu vs Re  |  %s  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', ...
              flow_model, Pr, Tb, Tw))
grid on

subplot(2,1,2)
plot(sweep.Re_bulk, sweep.Cf, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4)
xlabel('Re_{bulk}')
ylabel('C_f (Fanning) [-]')
title(sprintf('Cf vs Re  |  %s  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', ...
              flow_model, Pr, Tb, Tw))
grid on

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
plot(sweep_tur.Re_bulk, sweep_tur.Nu, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4)
xlabel('Re_{bulk}')
ylabel('Nusselt number [-]')
title(sprintf('Nu vs Re  |  turbulent  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', ...
    Pr, Tb, Tw))
grid on

subplot(2,1,2)
plot(sweep_tur.Re_bulk, sweep_tur.Cf, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4)
xlabel('Re_{bulk}')
ylabel('C_f (Fanning) [-]')
title(sprintf('Cf vs Re  |  turbulent  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', ...
    Pr, Tb, Tw))
grid on



%% ---- Plot combined regimes --- %%

figure('Name', 'thermoturb sweep — both regimes', 'NumberTitle', 'off');

subplot(2,1,1)
hold on
plot(sweep.Re_bulk,     sweep.Nu,     'b-o', 'LineWidth', 1.5, ...
     'MarkerSize', 4, 'DisplayName', 'Laminar')
plot(sweep_tur.Re_bulk, sweep_tur.Nu, 'r-o', 'LineWidth', 1.5, ...
     'MarkerSize', 4, 'DisplayName', 'Turbulent')
xlabel('Re_{bulk}')
ylabel('Nusselt number [-]')
title(sprintf('Nu vs Re  |  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', Pr, Tb, Tw))
legend('Location', 'best')
grid on

subplot(2,1,2)
hold on
plot(sweep.Re_bulk,     sweep.Cf,     'b-o', 'LineWidth', 1.5, ...
     'MarkerSize', 4, 'DisplayName', 'Laminar')
plot(sweep_tur.Re_bulk, sweep_tur.Cf, 'r-o', 'LineWidth', 1.5, ...
     'MarkerSize', 4, 'DisplayName', 'Turbulent')
xlabel('Re_{bulk}')
ylabel('C_f (Fanning) [-]')
title(sprintf('Cf vs Re  |  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', Pr, Tb, Tw))
legend('Location', 'best')
grid on