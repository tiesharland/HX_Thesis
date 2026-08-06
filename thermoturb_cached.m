function [Re_out, Nu_out, Cf_out] = thermoturb_cached(Re_bulk, Pr, Tb, Tw)

    persistent tt tt_loaded combos_lam combos_tur interp_built Re_tur_actual_min

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

    if isempty(interp_built) || ~interp_built

        for regime = [0, 1]
            sub  = tt(tt.regime == regime, :);
            grps = unique(sub{:, {'Pr','Tb','Tw'}}, 'rows');
            n    = size(grps, 1);

            entry.Pr = grps(:,1);
            entry.Tb = grps(:,2);
            entry.Tw = grps(:,3);
            entry.Re = cell(n,1);
            entry.Nu = cell(n,1);
            entry.Cf = cell(n,1);

            for i = 1:n
                mask = abs(sub.Pr - grps(i,1)) < 1e-4 & ...
                       abs(sub.Tb - grps(i,2)) < 0.1  & ...
                       abs(sub.Tw - grps(i,3)) < 0.1;
                g = sub(mask, :);
                [Re_sorted, idx] = sort(g.Re);
                entry.Re{i} = Re_sorted;
                entry.Nu{i} = g.Nu(idx);
                entry.Cf{i} = g.Cf(idx);
            end

            if regime == 0
                combos_lam = entry;
            else
                combos_tur = entry;
            end
            clear entry
        end

        % Actual minimum Re in turbulent data across all groups
        Re_tur_actual_min = min(cellfun(@min, combos_tur.Re));

        % Cap at Re_tur_min in case turbulent data extends below 3000
        Re_tur_actual_min = max(Re_tur_actual_min, Re_tur_min);

        interp_built = true;
        fprintf('thermoturb_cached: group structure built.\n');
        fprintf('  Laminar groups:         %d\n', numel(combos_lam.Pr));
        fprintf('  Turbulent groups:       %d\n', numel(combos_tur.Pr));
        fprintf('  Turbulent DNS Re start: %.0f\n', Re_tur_actual_min);
    end

    %% ---- Helper: lookup Nu at query (Re, Pr, Tb, Tw) ---- %%

    function val = lookup_Nu(Re_q, combos)
        n_grp        = numel(combos.Pr);
        Nu_at_groups = zeros(n_grp, 1);
        for i = 1:n_grp
            % No clamping — let interp1 extrapolate linearly beyond data range
            Nu_at_groups(i) = interp1(combos.Re{i}, combos.Nu{i}, Re_q, ...
                                      'linear', 'extrap');
        end
        Pr_q = max(min(combos.Pr), min(max(combos.Pr), Pr));
        Tb_q = max(min(combos.Tb), min(max(combos.Tb), Tb));
        Tw_q = max(min(combos.Tw), min(max(combos.Tw), Tw));
        F    = scatteredInterpolant(combos.Pr, combos.Tb, combos.Tw, Nu_at_groups, ...
                                    'linear', 'nearest');
        val  = F(Pr_q, Tb_q, Tw_q);
    end

    function val = lookup_Cf(Re_q, combos)
        n_grp        = numel(combos.Pr);
        Cf_at_groups = zeros(n_grp, 1);
        for i = 1:n_grp
            Cf_at_groups(i) = interp1(combos.Re{i}, combos.Cf{i}, Re_q, ...
                                       'linear', 'extrap');
        end
        Pr_q = max(min(combos.Pr), min(max(combos.Pr), Pr));
        Tb_q = max(min(combos.Tb), min(max(combos.Tb), Tb));
        Tw_q = max(min(combos.Tw), min(max(combos.Tw), Tw));
        F    = scatteredInterpolant(combos.Pr, combos.Tb, combos.Tw, Cf_at_groups, ...
                                    'linear', 'nearest');
        val  = F(Pr_q, Tb_q, Tw_q);
    end

    %% ---- Evaluate by regime ---- %%

    Re_out = Re_bulk;

    if Re_bulk <= Re_lam_max
        %--- Fully laminar: DNS laminar data ---%
        Nu_out = lookup_Nu(Re_bulk, combos_lam);
        Cf_out = lookup_Cf(Re_bulk, combos_lam);

    elseif Re_bulk >= Re_tur_min
        %--- Fully turbulent: DNS turbulent data ---%
        % Extrapolates linearly below Re_tur_actual_min down to Re_tur_min
        Nu_out = lookup_Nu(Re_bulk, combos_tur);
        Cf_out = lookup_Cf(Re_bulk, combos_tur);

    else
        %--- Transition 2300 < Re < 3000: cosine blend ---%
        % w = 0 at Re=2300 (fully laminar)
        % w = 1 at Re=3000 (fully turbulent)
        w = 0.5*(1 - cos(pi*(Re_bulk - Re_lam_max) / ...
                            (Re_tur_min - Re_lam_max)));

        Nu_lam = lookup_Nu(Re_bulk, combos_lam);
        Cf_lam = lookup_Cf(Re_bulk, combos_lam);
        Nu_tur = lookup_Nu(Re_bulk, combos_tur);
        Cf_tur = lookup_Cf(Re_bulk, combos_tur);

        Nu_out = (1-w)*Nu_lam + w*Nu_tur;
        Cf_out = (1-w)*Cf_lam + w*Cf_tur;

    end

    if isnan(Nu_out) || isnan(Cf_out)
        warning(['thermoturb_cached: NaN for Re=%.0f Pr=%.3f Tb=%.0f Tw=%.0f — ' ...
                 'inputs may be outside table coverage.'], Re_bulk, Pr, Tb, Tw);
    end

end