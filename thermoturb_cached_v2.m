function [Re_out, Nu_out, Cf_out] = thermoturb_cached(Re_bulk, Pr, Tb, Tw)

    persistent tt_loaded ...
               Re_lam_nodes F_Nu_lam F_Cf_lam ...
               Re_tur_nodes F_Nu_tur F_Cf_tur ...
               Re_tur_actual_min interp_built

    Re_lam_max = 2300;
    Re_tur_min = 3000;

    table_file = fullfile(pwd, 'thermoturb_table.mat');

    %% ---- Load and build if needed ---- %%

    file_info = dir(table_file);
    if isempty(tt_loaded) || file_info.datenum > tt_loaded || ...
       isempty(interp_built) || ~interp_built

        if ~isfile(table_file)
            error('thermoturb_cached: table not found. Run build_thermoturb_master_table.m first.');
        end

        loaded    = load(table_file, 'tt');
        tt        = loaded.tt;
        tt_loaded = file_info.datenum;

        fprintf('thermoturb_cached: building griddedInterpolant structure...\n');
        t_build = tic;

        for regime = [0, 1]
            sub  = tt(tt.regime == regime, :);

            Pr_u = unique(sub.Pr);
            Tb_u = unique(sub.Tb);
            Tw_u = unique(sub.Tw);
            nPr  = numel(Pr_u);
            nTb  = numel(Tb_u);
            nTw  = numel(Tw_u);

            grps  = unique(sub{:, {'Pr','Tb','Tw'}}, 'rows');
            n_grp = size(grps, 1);

            all_Re = unique(sub.Re);
            nRe    = numel(all_Re);

            Nu_grid = NaN(nRe, nPr, nTb, nTw);
            Cf_grid = NaN(nRe, nPr, nTb, nTw);

            if regime == 0
                regime_str = 'Laminar';
            else
                regime_str = 'Turbulent';
            end

            fprintf('  %s: filling %d groups x %d Re nodes...\n', ...
                    regime_str, n_grp, nRe);
            t_fill = tic;

            for i = 1:n_grp
                if mod(i, 100) == 0 || i == n_grp
                    fprintf('    Group %d/%d  (%.0f%%)  elapsed=%.1fs\n', ...
                            i, n_grp, 100*i/n_grp, toc(t_fill));
                end

                Pr_i = grps(i,1);
                Tb_i = grps(i,2);
                Tw_i = grps(i,3);

                mask = abs(sub.Pr - Pr_i) < 1e-4 & ...
                       abs(sub.Tb - Tb_i) < 0.1  & ...
                       abs(sub.Tw - Tw_i) < 0.1;
                g = sub(mask, :);
                [Re_sorted, idx] = sort(g.Re);

                iPr = find(abs(Pr_u - Pr_i) < 1e-4, 1);
                iTb = find(abs(Tb_u - Tb_i) < 0.5,  1);
                iTw = find(abs(Tw_u - Tw_i) < 0.5,  1);

                Nu_grid(:, iPr, iTb, iTw) = interp1(Re_sorted, g.Nu(idx), ...
                                                     all_Re, 'linear', 'extrap');
                Cf_grid(:, iPr, iTb, iTw) = interp1(Re_sorted, g.Cf(idx), ...
                                                     all_Re, 'linear', 'extrap');
            end

            fprintf('  %s grid filled in %.1fs.\n', regime_str, toc(t_fill));
            fprintf('  Building %d griddedInterpolants...\n', nRe);
            t_interp = tic;

            F_Nu = cell(nRe, 1);
            F_Cf = cell(nRe, 1);
            [Pr_grid_nd, Tb_grid_nd, Tw_grid_nd] = ndgrid(Pr_u, Tb_u, Tw_u);

            for iRe = 1:nRe
                if mod(iRe, 10) == 0 || iRe == nRe
                    fprintf('    Re node %d/%d  (%.0f%%)  elapsed=%.1fs\n', ...
                            iRe, nRe, 100*iRe/nRe, toc(t_interp));
                end
                F_Nu{iRe} = griddedInterpolant(Pr_grid_nd, Tb_grid_nd, Tw_grid_nd, ...
                                               squeeze(Nu_grid(iRe,:,:,:)), ...
                                               'linear', 'linear');
                F_Cf{iRe} = griddedInterpolant(Pr_grid_nd, Tb_grid_nd, Tw_grid_nd, ...
                                               squeeze(Cf_grid(iRe,:,:,:)), ...
                                               'linear', 'linear');
            end

            fprintf('  %s griddedInterpolants built in %.1fs.\n', ...
                    regime_str, toc(t_interp));

            if regime == 0
                Re_lam_nodes = all_Re;
                F_Nu_lam     = F_Nu;
                F_Cf_lam     = F_Cf;
            else
                Re_tur_nodes      = all_Re;
                F_Nu_tur          = F_Nu;
                F_Cf_tur          = F_Cf;
                Re_tur_actual_min = max(min(all_Re), Re_tur_min);
            end
        end

        interp_built = true;
        fprintf('thermoturb_cached: build complete in %.1fs.\n', toc(t_build));
        fprintf('  Laminar Re nodes:       %d  [%.0f - %.0f]\n', ...
                numel(Re_lam_nodes), min(Re_lam_nodes), max(Re_lam_nodes));
        fprintf('  Turbulent Re nodes:     %d  [%.0f - %.0f]\n', ...
                numel(Re_tur_nodes), min(Re_tur_nodes), max(Re_tur_nodes));
        fprintf('  Turbulent DNS Re start: %.0f\n', Re_tur_actual_min);
    end

    %% ---- Helper: lookup via pre-built griddedInterpolants + interp1 over Re ---- %%

    function [Nu_q, Cf_q] = lookup(Re_q, Re_nodes, F_Nu, F_Cf)
        nRe      = numel(Re_nodes);
        Nu_at_Re = zeros(nRe, 1);
        Cf_at_Re = zeros(nRe, 1);

        % Step 1: evaluate each Re-node griddedInterpolant at (Pr,Tb,Tw)
        % griddedInterpolant with 'linear','linear' extrapolates beyond grid
        for iRe = 1:nRe
            Nu_at_Re(iRe) = F_Nu{iRe}(Pr, Tb, Tw);
            Cf_at_Re(iRe) = F_Cf{iRe}(Pr, Tb, Tw);
        end

        % Step 2: linear interp/extrap over Re
        Nu_q = interp1(Re_nodes, Nu_at_Re, Re_q, 'linear', 'extrap');
        Cf_q = interp1(Re_nodes, Cf_at_Re, Re_q, 'linear', 'extrap');
    end

    %% ---- Evaluate by regime ---- %%

    Re_out = Re_bulk;

    if Re_bulk <= Re_lam_max
        %--- Fully laminar: DNS laminar data ---%
        [Nu_out, Cf_out] = lookup(Re_bulk, Re_lam_nodes, F_Nu_lam, F_Cf_lam);

    elseif Re_bulk >= Re_tur_min
        %--- Fully turbulent: DNS turbulent data ---%
        % Extrapolates linearly below Re_tur_actual_min down to Re_tur_min
        [Nu_out, Cf_out] = lookup(Re_bulk, Re_tur_nodes, F_Nu_tur, F_Cf_tur);

    else
        %--- Transition 2300 < Re < 3000: cosine blend ---%
        % w=0 at Re=2300 (fully laminar), w=1 at Re=3000 (fully turbulent)
        w = 0.5*(1 - cos(pi*(Re_bulk - Re_lam_max) / ...
                            (Re_tur_min - Re_lam_max)));

        [Nu_lam, Cf_lam] = lookup(Re_bulk, Re_lam_nodes, F_Nu_lam, F_Cf_lam);
        [Nu_tur, Cf_tur] = lookup(Re_bulk, Re_tur_nodes, F_Nu_tur, F_Cf_tur);

        Nu_out = (1-w)*Nu_lam + w*Nu_tur;
        Cf_out = (1-w)*Cf_lam + w*Cf_tur;

    end

    if isnan(Nu_out) || isnan(Cf_out)
        warning(['thermoturb_cached: NaN for Re=%.0f Pr=%.3f Tb=%.0f Tw=%.0f — ' ...
                 'inputs may be outside table coverage.'], Re_bulk, Pr, Tb, Tw);
    end

end