%% extend_thermoturb_table_lowTw.m
% Third extension of thermoturb_table.mat: extends the wall temperature
% range DOWN from 320 K to 300 K.
%
% New cells: Tw = 300:2:318, Tb = 260:2:Tw (so each new Tw row runs from
% Tb = 260 up to and including the diagonal Tb = Tw), at Pr 0.70 and 0.72.
% Tw >= 320 is already covered by the earlier extensions. The (Tb,Tw) hull
% of the table then has corners (260,300), (260,360), (360,360), (300,300).
%
% The existing table is loaded and only ever appended to. Pr values and Re
% ranges are the same as in the earlier extension scripts.
%
% IMPORTANT: do not run this at the same time as another extension script.
% Each script loads the table at the start and saves the whole table after
% every query, so two running together overwrite each other's new rows.
%
% Each (Pr,Tb,Tw,regime) is checked against the table only so a restart
% after an interruption doesn't append duplicate rows.

clear; clc;

%% ---- Definition of the new points ---- %%

Pr_values = [0.70 0.72];
Tw_values = 300 : 2 : 318;
Tb_min    = 260;

Re_lam_min =   10;   Re_lam_max =  5000;
Re_tur_min = 5400;   Re_tur_max = 10000;

% Only used for the restart check
tol_Pr = 0.005;
tol_Tb = 1;
tol_Tw = 1;

%% ---- Build the explicit list of cells to query ---- %%

cells = zeros(0, 3);   % [Pr Tb Tw]
for Pr = Pr_values
    for Tw = Tw_values
        for Tb = Tb_min : 2 : Tw
            cells(end+1, :) = [Pr Tb Tw]; %#ok<SAGROW>
        end
    end
end
n_cells = size(cells, 1);
n_total = 2 * n_cells;   % laminar + turbulent each
fprintf('%d (Pr,Tb,Tw) cells x 2 regimes = %d queries (~%.1f h at 6 s each)\n', ...
        n_cells, n_total, n_total*6/3600);

%% ---- Load existing master table (never replaced, only appended to) ---- %%

table_file = fullfile(pwd, 'thermoturb_table.mat');
if ~isfile(table_file)
    error('extend_thermoturb_table: %s not found -- run from the folder holding the existing table.', table_file);
end
loaded = load(table_file, 'tt');
tt     = loaded.tt;
fprintf('Loaded existing table with %d entries.\n', height(tt));

%% ---- Query loop ---- %%

regime_names = {'laminar', 'turbulent'};
n_queried = 0; n_skipped = 0; n_failed = 0;
t_per_query = zeros(n_total, 1);
pause_s = 1.5;   % politeness pause between queries to thermoturb.com
t_start = tic;

for c = 1:n_cells
    Pr = cells(c,1); Tb = cells(c,2); Tw = cells(c,3);

    for regime_flag = 0:1
        flow_model = regime_names{regime_flag + 1};

        if regime_flag == 0
            Re_min = Re_lam_min; Re_max = Re_lam_max;
        else
            Re_min = Re_tur_min; Re_max = Re_tur_max;
        end

        % restart check only
        already = sum(abs(tt.Pr - Pr) < tol_Pr & abs(tt.Tb - Tb) < tol_Tb & ...
                      abs(tt.Tw - Tw) < tol_Tw & tt.regime == regime_flag) >= 5;
        if already
            n_skipped = n_skipped + 1;
            continue
        end

        fprintf('  -> Pr=%.2f Tb=%.0fK Tw=%.0fK %s Re=[%d-%d] ... ', ...
                Pr, Tb, Tw, flow_model, Re_min, Re_max);

        t_query = tic;
        try
            sweep = thermoturb_query_sweep(Re_min, Re_max, Pr, Tb, Tw, flow_model);

            t_per_query(n_queried+1) = toc(t_query);
            n_queried = n_queried + 1;

            n_new = height(sweep);
            new_rows = table( ...
                repmat(Pr, n_new, 1), repmat(Tb, n_new, 1), repmat(Tw, n_new, 1), ...
                sweep.Re_bulk, sweep.Nu, sweep.Cf, sweep.Stanton, ...
                sweep.Retau, sweep.Retau_cp, repmat(double(regime_flag), n_new, 1), ...
                'VariableNames', {'Pr','Tb','Tw','Re','Nu','Cf', ...
                                  'Stanton','Retau','Retau_cp','regime'});

            tt = [tt; new_rows]; %#ok<AGROW>
            save(table_file, 'tt');

            fprintf('OK -- %d points in %.1fs. Table: %d entries.\n', ...
                    n_new, t_per_query(n_queried), height(tt));
        catch ME
            n_failed = n_failed + 1;
            fprintf('FAILED in %.1fs: %s\n', toc(t_query), ME.message);
        end

        pause(pause_s);

        print_progress(n_queried + n_skipped + n_failed, n_total, n_failed, ...
                       t_start, t_per_query, n_queried, pause_s);
    end
end

%% ---- Summary ---- %%

fprintf('\n========================================\n');
fprintf('Done in %.1f minutes.\n', toc(t_start)/60);
fprintf('  Queried: %d   Already present: %d   Failed: %d\n', n_queried, n_skipped, n_failed);
fprintf('  Table entries: %d\n', height(tt));
fprintf('========================================\n');
fprintf('Failed queries are not retried automatically -- re-run this script to retry them.\n');

%% ---- Local helper: one-line progress / ETA after each query ---- %%
% The ETA assumes every remaining query is actually sent (cells already in
% the table are skipped instantly, so after a restart it overestimates), and
% uses the average of the last 50 queries plus the pause, so it follows
% any change in server speed.
function print_progress(n_done, n_total, n_failed, t_start, t_per_query, n_queried, pause_s)
elapsed = toc(t_start);
if n_queried == 0
    fprintf('     progress %d/%d | %d failed | elapsed %s | ETA: no successful query yet\n', ...
            n_done, n_total, n_failed, fmt_duration(elapsed));
    return
end
recent  = t_per_query(max(1, n_queried - 49):n_queried);
eta_s   = (n_total - n_done) * (mean(recent) + pause_s);
finish  = datetime('now') + seconds(eta_s);
finish.Format = 'eee dd-MMM HH:mm';
fprintf('     progress %d/%d (%.1f%%) | %d failed | elapsed %s | ETA %s | finish ~ %s\n', ...
        n_done, n_total, 100*n_done/n_total, n_failed, ...
        fmt_duration(elapsed), fmt_duration(eta_s), char(finish));
end

function s = fmt_duration(sec)
if sec >= 3600
    s = sprintf('%dh %02dm', floor(sec/3600), floor(mod(sec, 3600)/60));
else
    s = sprintf('%dm %02ds', floor(sec/60), floor(mod(sec, 60)));
end
end