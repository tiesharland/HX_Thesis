%% build_thermoturb_master_table.m

clear; clc;

%% ---- Grid definition ---- %%

Pr_values = 0.700 : 0.002 : 0.720;   % 11 values
Tb_values = 264   : 4     : 330;      % 17 valuesuni
Tw_values = 340   : 2     : 358;      % 10 values
% Total: 11 x 17 x 10 x 2 regimes = 3740 queries
% At ~4.5s per query -> ~4.7 hours

Re_lam_min =   10;   Re_lam_max =  5000;
Re_tur_min = 5400;   Re_tur_max = 10000;

tol_Pr = 0.001;
tol_Tb = 2;     % K — half the step size
tol_Tw = 1;     % K — half the step size

%% ---- Load or initialise master table ---- %%

table_file = fullfile(pwd, 'thermoturb_table.mat');

if isfile(table_file)
    loaded = load(table_file, 'tt');
    tt     = loaded.tt;
    fprintf('Loaded existing table with %d entries.\n', height(tt));
else
    tt = table('Size', [0 10], ...
               'VariableTypes', repmat({'double'}, 1, 10), ...
               'VariableNames', {'Pr','Tb','Tw','Re','Nu','Cf', ...
                                 'Stanton','Retau','Retau_cp','regime'});
    save(table_file, 'tt');
    fprintf('Initialised new empty master table.\n');
end

%% ---- Summary ---- %%

n_Pr    = numel(Pr_values);
n_Tb    = numel(Tb_values);
n_Tw    = numel(Tw_values);
n_total = n_Pr * n_Tb * n_Tw * 2;
fprintf('\nGrid: %d Pr x %d Tb x %d Tw x 2 regimes = %d queries\n', ...
        n_Pr, n_Tb, n_Tw, n_total);
fprintf('Estimated time: %.1f hours (assuming 4.5s/query)\n\n', ...
        n_total*4.5/3600);

%% ---- Batch query loop ---- %%

n_queried  = 0;
n_skipped  = 0;
n_failed   = 0;
query_num  = 0;
t_start    = tic;
t_per_query = zeros(n_total, 1);

for i_Pr = 1:n_Pr
for i_Tb = 1:n_Tb
for i_Tw = 1:n_Tw

    Pr = Pr_values(i_Pr);
    Tb = Tb_values(i_Tb);
    Tw = Tw_values(i_Tw);

    if Tw <= Tb
        n_skipped = n_skipped + 2;
        continue
    end

    if height(tt) > 0
        mask_base = abs(tt.Pr - Pr) < tol_Pr & ...
                    abs(tt.Tb - Tb) < tol_Tb  & ...
                    abs(tt.Tw - Tw) < tol_Tw;
        has_lam   = sum(mask_base & tt.regime == 0) >= 5;
        has_tur   = sum(mask_base & tt.regime == 1) >= 5;
    else
        has_lam = false;
        has_tur = false;
    end

    regimes_to_query = {};
    if ~has_lam, regimes_to_query{end+1} = 'laminar';   end %#ok<AGROW>
    if ~has_tur, regimes_to_query{end+1} = 'turbulent'; end %#ok<AGROW>

    if isempty(regimes_to_query)
        n_skipped = n_skipped + 2;
        continue
    end

    for i_regime = 1:numel(regimes_to_query)

        flow_model = regimes_to_query{i_regime};
        query_num  = query_num + 1;

        % n_done_total counts both queried and skipped so far
        n_done_total    = n_queried + n_skipped + n_failed;
        n_remaining_tot = n_total - n_done_total;

        if strcmp(flow_model, 'laminar')
            Re_min      = Re_lam_min;
            Re_max      = Re_lam_max;
            regime_flag = 0;
        else
            Re_min      = Re_tur_min;
            Re_max      = Re_tur_max;
            regime_flag = 1;
        end

        % Time estimate based on average of actual queries only
        % (skipped ones take negligible time so excluded from avg)
        if n_queried > 0
            avg_t     = mean(t_per_query(1:n_queried));
            % Only unqueried remaining combinations need time
            n_still_to_query = n_total - n_done_total - ...
                               (numel(regimes_to_query) - i_regime);
            remaining = max(0, n_still_to_query) * avg_t;
            fprintf('[%d/%d queried | %d/%d skipped | %d remaining] ', ...
                    n_queried, n_total-n_skipped, n_skipped, n_total, n_remaining_tot);
            fprintf('avg=%.1fs  elapsed=%.0fm  remaining~%.0fm\n', ...
                    avg_t, toc(t_start)/60, remaining/60);
        else
            fprintf('[%d/%d queried | %d/%d skipped | %d remaining] ', ...
                    n_queried, n_total-n_skipped, n_skipped, n_total, n_remaining_tot);
            fprintf('elapsed=%.0fm  remaining: estimating...\n', toc(t_start)/60);
        end

        fprintf('  -> Pr=%.3f Tb=%.0fK Tw=%.0fK %s Re=[%d-%d] ... ', ...
                Pr, Tb, Tw, flow_model, Re_min, Re_max);

        t_query = tic;

        try
            sweep = thermoturb_query_sweep(Re_min, Re_max, Pr, Tb, Tw, flow_model);

            t_per_query(n_queried+1) = toc(t_query);
            n_queried = n_queried + 1;

            n_new    = height(sweep);
            new_rows = table(...
                repmat(Pr,                   n_new, 1), ...
                repmat(Tb,                   n_new, 1), ...
                repmat(Tw,                   n_new, 1), ...
                sweep.Re_bulk,                          ...
                sweep.Nu,                               ...
                sweep.Cf,                               ...
                sweep.Stanton,                          ...
                sweep.Retau,                            ...
                sweep.Retau_cp,                         ...
                repmat(double(regime_flag),  n_new, 1), ...
                'VariableNames', {'Pr','Tb','Tw','Re','Nu','Cf', ...
                                  'Stanton','Retau','Retau_cp','regime'});

            tt = [tt; new_rows]; %#ok<AGROW>
            save(table_file, 'tt');

            fprintf('OK — %d points in %.1fs. Table: %d entries.\n', ...
                    n_new, t_per_query(n_queried), height(tt));

        catch ME
            t_per_query(n_queried+1) = toc(t_query);
            fprintf('FAILED in %.1fs: %s\n', t_per_query(n_queried+1), ME.message);
            n_failed = n_failed + 1;
        end

        pause(1.5);

    end

end
end
end

%% ---- Final summary ---- %%

total_time = toc(t_start);
fprintf('\n========================================\n');
fprintf('Build complete in %.1f minutes.\n', total_time/60);
fprintf('  Queried:       %d\n', n_queried);
fprintf('  Skipped:       %d\n', n_skipped);
fprintf('  Failed:        %d\n', n_failed);
fprintf('  Table entries: %d\n', height(tt));
if n_queried > 0
    fprintf('  Avg time/query: %.2fs\n', mean(t_per_query(1:n_queried)));
end
fprintf('  Saved to:      %s\n', table_file);
fprintf('========================================\n');
fprintf('\nTo refine later, reduce Tb_step to 2 and re-run.\n');
fprintf('Existing entries will be skipped automatically.\n');