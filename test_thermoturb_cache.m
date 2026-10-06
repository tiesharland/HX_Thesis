%% test_thermoturb_cache.m
% Plots Nu and Cf from thermoturb_cached over the full Re range,
% showing laminar, transition, extrapolated turbulent, and DNS turbulent
% regimes continuously. Infers turbulent data boundary from table directly.

clear; clc;

%% ---- Fixed conditions ---- %%

Pr = 0.703;
Tb = 330;    % K — bulk air temperature
Tw = 356;    % K — wall temperature

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

%% ---- Re sweep ---- %%

Re_lam_max  = 2300;
Re_tur_min  = 3000;
Re_max_plot = 10000;

Re_arr = unique([...
    linspace(100,               Re_lam_max,        50), ...
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
    [~, Nu_arr(i), Cf_arr(i)] = thermoturb_cached(Re_arr(i), Pr, Tb, Tw);
end
fprintf('Done in %.2fs.\n\n', toc);

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

%% ---- Colour scheme ---- %%

col_lam = [0.122 0.471 0.706];   % blue   — laminar DNS
col_tra = [0.839 0.153 0.157];   % red    — transition cosine blend
col_ext = [0.580 0.404 0.741];   % purple — turbulent extrapolated
col_tur = [0.173 0.627 0.173];   % green  — turbulent DNS
lw = 2;

%% ---- Figure 1: Nu vs Re ---- %%

figure('Name', 'Nu vs Re — thermoturb_cached', 'NumberTitle', 'off');
hold on

plot(Re_arr(is_lam), Nu_arr(is_lam), '-',  'Color', col_lam, 'LineWidth', lw, ...
     'DisplayName', 'Laminar (DNS)')
plot(Re_arr(is_tra), Nu_arr(is_tra), '--', 'Color', col_tra, 'LineWidth', lw, ...
     'DisplayName', 'Transition (cosine blend)')
if any(is_ext)
    plot(Re_arr(is_ext), Nu_arr(is_ext), ':', 'Color', col_ext, 'LineWidth', lw, ...
         'DisplayName', sprintf('Turbulent (extrap. from Re=%.0f)', Re_tur_actual_min))
end
plot(Re_arr(is_tur), Nu_arr(is_tur), '-',  'Color', col_tur, 'LineWidth', lw, ...
     'DisplayName', sprintf('Turbulent DNS (Re\\geq%.0f)', Re_tur_actual_min))

xline(Re_lam_max, '--k', 'Re=2300', 'LabelHorizontalAlignment', 'left',  'LineWidth', 1)
xline(Re_tur_min, '--k', 'Re=3000', 'LabelHorizontalAlignment', 'right', 'LineWidth', 1)
if Re_tur_actual_min > Re_tur_min
    xline(Re_tur_actual_min, '--k', sprintf('Re=%.0f (data)', Re_tur_actual_min), ...
          'LabelHorizontalAlignment', 'right', 'LineWidth', 1)
end

xlabel('Re_{bulk}')
ylabel('Nusselt number [-]')
title(sprintf('Nu vs Re  |  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', Pr, Tb, Tw))
legend('Location', 'northwest')
grid on

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

xline(Re_lam_max, '--k', 'Re=2300', 'LabelHorizontalAlignment', 'left',  'LineWidth', 1)
xline(Re_tur_min, '--k', 'Re=3000', 'LabelHorizontalAlignment', 'right', 'LineWidth', 1)
if Re_tur_actual_min > Re_tur_min
    xline(Re_tur_actual_min, '--k', sprintf('Re=%.0f (data)', Re_tur_actual_min), ...
          'LabelHorizontalAlignment', 'right', 'LineWidth', 1)
end

xlabel('Re_{bulk}')
ylabel('C_f (Fanning) [-]')
title(sprintf('Cf vs Re  |  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', Pr, Tb, Tw))
legend('Location', 'northeast')
grid on

%% ---- Figure 3: Combined 2x1 panel ---- %%

figure('Name', 'Nu and Cf vs Re — thermoturb_cached', 'NumberTitle', 'off');

subplot(2,1,1)
hold on
plot(Re_arr(is_lam), Nu_arr(is_lam), '-',  'Color', col_lam, 'LineWidth', lw, ...
     'DisplayName', 'Laminar (DNS)')
plot(Re_arr(is_tra), Nu_arr(is_tra), '--', 'Color', col_tra, 'LineWidth', lw, ...
     'DisplayName', 'Transition')
if any(is_ext)
    plot(Re_arr(is_ext), Nu_arr(is_ext), ':', 'Color', col_ext, 'LineWidth', lw, ...
         'DisplayName', 'Turbulent (extrap.)')
end
plot(Re_arr(is_tur), Nu_arr(is_tur), '-',  'Color', col_tur, 'LineWidth', lw, ...
     'DisplayName', 'Turbulent (DNS)')
xline(Re_lam_max, '--k', 'Re=2300', 'LabelHorizontalAlignment', 'left',  'LineWidth', 1, 'HandleVisibility', 'off')
xline(Re_tur_min, '--k', 'Re=3000', 'LabelHorizontalAlignment', 'right', 'LineWidth', 1, 'HandleVisibility', 'off')
if Re_tur_actual_min > Re_tur_min
    xline(Re_tur_actual_min, '--k', sprintf('Re=%.0f', Re_tur_actual_min), ...
          'LabelHorizontalAlignment', 'right', 'LineWidth', 1, 'HandleVisibility', 'off')
end
ylabel('Nusselt number [-]')
title(sprintf('Nu and Cf vs Re  |  Pr=%.3f  Tb=%.0fK  Tw=%.0fK', Pr, Tb, Tw))
legend('Location', 'northwest', 'FontSize', 8)
grid on

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
xlabel('Re_{bulk}')
ylabel('C_f (Fanning) [-]')
legend('Location', 'northeast', 'FontSize', 8)
grid on

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