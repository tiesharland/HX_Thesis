function [Re_out, Nu_out, Cf_out] = thermoturb_cached(Re_bulk, Pr, Tb, Tw)

    % Set true to log every query (Re/Pr/Tb/Tw, regime, and whether each
    % dimension fell outside the DNS table's coverage) via
    % thermoturb_query_log, for later coverage/extrapolation plotting with
    % plot_thermoturb_coverage.m. Call thermoturb_query_log('reset') before
    % a run if you only want that run's queries in the log.
    log_queries = true;

    persistent tt tt_loaded interp_built Re_tur_actual_min ...
               Re_lam_rng Re_tur_rng Pr_rng_lam Tb_rng_lam Tw_rng_lam ...
               Pr_rng_tur Tb_rng_tur Tw_rng_tur ...
               Pr_lam Tb_lam Tw_lam Pr_tur Tb_tur Tw_tur ...
               Re_mat_lam Nu_mat_lam Cf_mat_lam npts_lam ...
               Re_mat_tur Nu_mat_tur Cf_mat_tur npts_tur ...
               F_Nu_lam F_Cf_lam F_Nu_tur F_Cf_tur

    Re_lam_max = 2300;
    Re_tur_min = 3000;

    table_file = fullfile(pwd, 'thermoturb_table.mat');

    %% ---- Load table if not yet loaded or file has changed ---- %%

    file_info = dir(table_file);
    if isempty(tt_loaded) || file_info.datenum > tt_loaded
        if ~isfile(table_file)
            error('thermoturb_cached: table not found. Run build_thermoturb_master_table.m first.');
        end
        loaded       = load(table_file, 'tt');
        tt           = loaded.tt;
        tt_loaded    = file_info.datenum;
        interp_built = false;
        fprintf('thermoturb_cached: table loaded — %d entries.\n', height(tt));
    end

    %% ---- Build persistent group structure if needed ---- %%
    %
    % Groups are stored as dense, Re-sorted, Inf/NaN-padded matrices
    % (Re_mat/Nu_mat/Cf_mat, n_groups x max_points_per_group) instead of
    % per-group cells, so every query's Re-interpolation step runs as one
    % vectorized pass over all groups at once rather than a MATLAB for-loop
    % calling interp1 once per group (there can be ~1000+ groups, and this
    % runs on every thermoturb_cached call, so the loop overhead matters).
    % Re is padded with +Inf (never counts as "<= Re_q") and Nu/Cf are
    % padded with NaN (never read, since the padded columns are always
    % beyond each row's clamped interpolation index).

    if isempty(interp_built) || ~interp_built

        for regime = [0, 1]
            sub  = tt(tt.regime == regime, :);
            grps = unique(sub{:, {'Pr','Tb','Tw'}}, 'rows');
            n    = size(grps, 1);

            Pr_g = grps(:,1); Tb_g = grps(:,2); Tw_g = grps(:,3);
            Re_cell = cell(n,1); Nu_cell = cell(n,1); Cf_cell = cell(n,1);
            npts    = zeros(n,1);

            for i = 1:n
                mask = abs(sub.Pr - grps(i,1)) < 1e-4 & ...
                       abs(sub.Tb - grps(i,2)) < 0.1  & ...
                       abs(sub.Tw - grps(i,3)) < 0.1;
                g = sub(mask, :);
                [Re_sorted, idx] = sort(g.Re);
                Re_cell{i} = Re_sorted;
                Nu_cell{i} = g.Nu(idx);
                Cf_cell{i} = g.Cf(idx);
                npts(i)    = numel(Re_sorted);
            end

            Mmax   = max(npts);
            Re_mat = inf(n, Mmax);
            Nu_mat = nan(n, Mmax);
            Cf_mat = nan(n, Mmax);
            for i = 1:n
                Re_mat(i, 1:npts(i)) = Re_cell{i};
                Nu_mat(i, 1:npts(i)) = Nu_cell{i};
                Cf_mat(i, 1:npts(i)) = Cf_cell{i};
            end

            if regime == 0
                Pr_lam = Pr_g; Tb_lam = Tb_g; Tw_lam = Tw_g;
                Re_mat_lam = Re_mat; Nu_mat_lam = Nu_mat; Cf_mat_lam = Cf_mat;
                npts_lam   = npts;
            else
                Pr_tur = Pr_g; Tb_tur = Tb_g; Tw_tur = Tw_g;
                Re_mat_tur = Re_mat; Nu_mat_tur = Nu_mat; Cf_mat_tur = Cf_mat;
                npts_tur   = npts;
            end
        end

        % Actual minimum Re in turbulent data across all groups (column 1
        % is always each row's minimum, since rows are Re-sorted ascending
        % and never padded at the start).
        Re_tur_actual_min = min(Re_mat_tur(:,1));

        % Cap at Re_tur_min in case turbulent data extends below 3000
        Re_tur_actual_min = max(Re_tur_actual_min, Re_tur_min);

        % Coverage bounds per regime, for flagging extrapolation/clamping
        % of logged queries (used by plot_thermoturb_coverage.m). Each
        % row's max is at column npts(i) (last valid, pre-padding entry).
        n_lam        = numel(Pr_lam);
        n_tur        = numel(Pr_tur);
        lam_last_idx = (1:n_lam)' + (npts_lam - 1) * n_lam;
        tur_last_idx = (1:n_tur)' + (npts_tur - 1) * n_tur;
        Re_lam_rng   = [min(Re_mat_lam(:,1)), max(Re_mat_lam(lam_last_idx))];
        Re_tur_rng   = [min(Re_mat_tur(:,1)), max(Re_mat_tur(tur_last_idx))];
        Pr_rng_lam   = [min(Pr_lam), max(Pr_lam)];
        Tb_rng_lam   = [min(Tb_lam), max(Tb_lam)];
        Tw_rng_lam   = [min(Tw_lam), max(Tw_lam)];
        Pr_rng_tur   = [min(Pr_tur), max(Pr_tur)];
        Tb_rng_tur   = [min(Tb_tur), max(Tb_tur)];
        Tw_rng_tur   = [min(Tw_tur), max(Tw_tur)];

        % Build the (Pr,Tb,Tw) triangulation ONCE per regime. The
        % triangulation never changes — only the Nu/Cf values at these
        % fixed points change per query (they're Re-interpolated first).
        % scatteredInterpolant lets us swap .Values on an existing object
        % without re-triangulating, which is the expensive part — doing
        % that from scratch on every single query (the old behavior) is
        % what was making long discretised/DNS runs slow.
        F_Nu_lam = scatteredInterpolant(Pr_lam, Tb_lam, Tw_lam, zeros(n_lam,1), 'linear', 'nearest');
        F_Cf_lam = scatteredInterpolant(Pr_lam, Tb_lam, Tw_lam, zeros(n_lam,1), 'linear', 'nearest');
        F_Nu_tur = scatteredInterpolant(Pr_tur, Tb_tur, Tw_tur, zeros(n_tur,1), 'linear', 'nearest');
        F_Cf_tur = scatteredInterpolant(Pr_tur, Tb_tur, Tw_tur, zeros(n_tur,1), 'linear', 'nearest');

        interp_built = true;
        fprintf('thermoturb_cached: group structure built.\n');
        fprintf('  Laminar groups:         %d\n', n_lam);
        fprintf('  Turbulent groups:       %d\n', n_tur);
        fprintf('  Turbulent DNS Re start: %.0f\n', Re_tur_actual_min);
    end

    %% ---- Helper: vectorized Re-interpolation across all groups at once ---- %%
    % Re_mat is Inf-padded past each row's npts(i) valid columns, so
    % summing (Re_mat <= Re_q) along rows gives, for each group, how many
    % of its real sample points sit at or below the query Re — exactly the
    % bracketing index linear interpolation needs. Clamping that count to
    % [1, npts-1] reproduces interp1(...,'linear','extrap')'s behavior:
    % inside the group's range it picks the right bracketing segment, and
    % outside it, it keeps extending the nearest edge segment's slope.
    function [Nu_vals, Cf_vals] = interp_groups(Re_q, Re_mat, Nu_mat, Cf_mat, npts)
        n_grp = size(Re_mat, 1);
        idx   = sum(Re_mat <= Re_q, 2);
        idx   = max(1, min(idx, npts - 1));
        rows  = (1:n_grp)';
        lin0  = rows + (idx - 1) * n_grp;
        lin1  = rows + idx * n_grp;
        x0    = Re_mat(lin0);
        x1    = Re_mat(lin1);
        t     = (Re_q - x0) ./ (x1 - x0);
        Nu_vals = Nu_mat(lin0) + t .* (Nu_mat(lin1) - Nu_mat(lin0));
        Cf_vals = Cf_mat(lin0) + t .* (Cf_mat(lin1) - Cf_mat(lin0));
    end

    %% ---- Helper: full lookup (Re-interp across groups + Pr/Tb/Tw interp) ---- %%
    % F_Nu/F_Cf are pre-triangulated scatteredInterpolants for this regime
    % (built once above); only their .Values get updated here, so no
    % retriangulation happens on the hot path.
    function [Nu_val, Cf_val] = lookup_pair(Re_q, Re_mat, Nu_mat, Cf_mat, npts, ...
                                             Pr_rng, Tb_rng, Tw_rng, F_Nu, F_Cf)
        [Nu_at_groups, Cf_at_groups] = interp_groups(Re_q, Re_mat, Nu_mat, Cf_mat, npts);
        Pr_q = max(Pr_rng(1), min(Pr_rng(2), Pr));
        Tb_q = max(Tb_rng(1), min(Tb_rng(2), Tb));
        Tw_q = max(Tw_rng(1), min(Tw_rng(2), Tw));
        F_Nu.Values = Nu_at_groups;
        F_Cf.Values = Cf_at_groups;
        Nu_val = F_Nu(Pr_q, Tb_q, Tw_q);
        Cf_val = F_Cf(Pr_q, Tb_q, Tw_q);
    end

    %% ---- Evaluate by regime ---- %%

    Re_out = Re_bulk;

    if Re_bulk <= Re_lam_max
        %--- Fully laminar: DNS laminar data ---%
        [Nu_out, Cf_out] = lookup_pair(Re_bulk, Re_mat_lam, Nu_mat_lam, Cf_mat_lam, npts_lam, ...
                                        Pr_rng_lam, Tb_rng_lam, Tw_rng_lam, F_Nu_lam, F_Cf_lam);

        regime_tag = 0;
        Re_extrap  = Re_bulk < Re_lam_rng(1) || Re_bulk > Re_lam_rng(2);
        Pr_extrap  = Pr < Pr_rng_lam(1) || Pr > Pr_rng_lam(2);
        Tb_extrap  = Tb < Tb_rng_lam(1) || Tb > Tb_rng_lam(2);
        Tw_extrap  = Tw < Tw_rng_lam(1) || Tw > Tw_rng_lam(2);

    elseif Re_bulk >= Re_tur_min
        %--- Fully turbulent: DNS turbulent data ---%
        % Extrapolates linearly below Re_tur_actual_min down to Re_tur_min
        [Nu_out, Cf_out] = lookup_pair(Re_bulk, Re_mat_tur, Nu_mat_tur, Cf_mat_tur, npts_tur, ...
                                        Pr_rng_tur, Tb_rng_tur, Tw_rng_tur, F_Nu_tur, F_Cf_tur);

        regime_tag = 1;
        Re_extrap  = Re_bulk < Re_tur_rng(1) || Re_bulk > Re_tur_rng(2);
        Pr_extrap  = Pr < Pr_rng_tur(1) || Pr > Pr_rng_tur(2);
        Tb_extrap  = Tb < Tb_rng_tur(1) || Tb > Tb_rng_tur(2);
        Tw_extrap  = Tw < Tw_rng_tur(1) || Tw > Tw_rng_tur(2);

    else
        %--- Transition 2300 < Re < 3000: cosine blend ---%
        % w = 0 at Re=2300 (fully laminar)
        % w = 1 at Re=3000 (fully turbulent)
        w = 0.5*(1 - cos(pi*(Re_bulk - Re_lam_max) / ...
                            (Re_tur_min - Re_lam_max)));

        [Nu_lam, Cf_lam] = lookup_pair(Re_bulk, Re_mat_lam, Nu_mat_lam, Cf_mat_lam, npts_lam, ...
                                        Pr_rng_lam, Tb_rng_lam, Tw_rng_lam, F_Nu_lam, F_Cf_lam);
        [Nu_tur, Cf_tur] = lookup_pair(Re_bulk, Re_mat_tur, Nu_mat_tur, Cf_mat_tur, npts_tur, ...
                                        Pr_rng_tur, Tb_rng_tur, Tw_rng_tur, F_Nu_tur, F_Cf_tur);

        Nu_out = (1-w)*Nu_lam + w*Nu_tur;
        Cf_out = (1-w)*Cf_lam + w*Cf_tur;

        % Re itself is by definition inside the designed blend zone, so it
        % isn't flagged as extrapolated here; Pr/Tb/Tw are flagged if
        % outside EITHER regime's coverage box, since both were evaluated.
        regime_tag = 2;
        Re_extrap  = false;
        Pr_extrap  = Pr < Pr_rng_lam(1) || Pr > Pr_rng_lam(2) || ...
                     Pr < Pr_rng_tur(1) || Pr > Pr_rng_tur(2);
        Tb_extrap  = Tb < Tb_rng_lam(1) || Tb > Tb_rng_lam(2) || ...
                     Tb < Tb_rng_tur(1) || Tb > Tb_rng_tur(2);
        Tw_extrap  = Tw < Tw_rng_lam(1) || Tw > Tw_rng_lam(2) || ...
                     Tw < Tw_rng_tur(1) || Tw > Tw_rng_tur(2);

    end

    if log_queries
        entry.Re        = Re_bulk;
        entry.Pr        = Pr;
        entry.Tb        = Tb;
        entry.Tw        = Tw;
        entry.regime    = regime_tag;
        entry.Re_extrap = Re_extrap;
        entry.Pr_extrap = Pr_extrap;
        entry.Tb_extrap = Tb_extrap;
        entry.Tw_extrap = Tw_extrap;
        thermoturb_query_log('log', entry);
    end

    if isnan(Nu_out) || isnan(Cf_out)
        warning(['thermoturb_cached: NaN for Re=%.0f Pr=%.3f Tb=%.0f Tw=%.0f — ' ...
                 'inputs may be outside table coverage.'], Re_bulk, Pr, Tb, Tw);
    end

end