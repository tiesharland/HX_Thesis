function comparison = compare_cfd_meshes_vs_2d(cfd, varargin)
% compare_cfd_meshes_vs_2d  Plot several CFD unit-cell results (e.g. the same
% case on a coarse / medium / fine / finer mesh, each from
% process_wall_CFD_data) against the 2D discretised model - with and/or
% without DNS - and the 1D discretised model (N_cool_seg = 0) on shared axes, one figure per quantity (air bulk
% temperature, static pressure, h, velocity, mass flow rate, Nu, and density
% when the CFD structs carry results.rho).
%
% Same plots and conventions as compare_cfd_vs_2d (each series at its own
% native x positions, no interpolation; the model is compared at coolant
% segment l = 1 by default; h/Nu leave out the z = 0 entrance station).
%
% Usage:
%   % 2D model results read from the results folder (whatever exists):
%   comparison = compare_cfd_meshes_vs_2d({cfd_coarse, cfd_medium, cfd_fine, cfd_finer}, ...
%       'cfd_labels', {'coarse','medium','fine','finer'}, ...
%       'id', 1, 'N_segments', 50, 'N_cool_seg', 5);
%
%   % or hand over the 2D results yourself:
%   comparison = compare_cfd_meshes_vs_2d(cfd_list, 'results_2D', {res_corr, res_dns}, ...
%       'labels_2D', {'2D correlations','2D DNS'});
%
% cfd           : cell array of structs from process_wall_CFD_data (a single
%                 struct also works). Give them from coarsest to finest: the
%                 CFD lines go from light to dark grey in that order.
% Name-value options:
%   'cfd_labels' : cell array of legend names for the CFD cases
%                  (default 'CFD 1', 'CFD 2', ...).
%   'results_2D' / 'labels_2D' : 2D results structs (and their legend names)
%                  supplied directly. Leave out to read from the results folder.
%   'id', 'N_segments', 'N_cool_seg' : which saved 2D case to read from
%                  results/id-<id>/ (N_cool_seg >= 1 for the 2D path).
%   'use_DNS'    : which variants to look for, default [false true]. A variant
%                  with no saved file is skipped with a warning; if none is
%                  found the function stops with an error. If several saved
%                  versions of a file exist (-1, -2, ...), the base file is used.
%   'include_1d' : also look for the 1D discretised model (N_cool_seg = 0, same
%                  id and N_segments) for every use_DNS variant (default true;
%                  missing ones are skipped with a warning). Not needed when
%                  results_2D is given -- put 1D results in that list instead.
%   'include_lumped' : also plot the lumped model (default true) if
%                  results/id-<id>/lumped_<id>.mat exists; it has no position
%                  resolution, so its values are drawn as straight lines from
%                  the HX inlet to the outlet (T and P: inlet -> outlet values;
%                  h, v, Nu: constant). No mass-flow curve for it.
%   'cfd_P_operating' : operating (reference) pressure of the CFD case [Pa], e.g. 58200.
%                  Fluent's Pbulk is gauge pressure, so the CFD static pressure is drawn
%                  as Pbulk + this value, i.e. the actual pressure, comparable to the
%                  model's P_air_seg. Leave empty to fall back to the old behaviour
%                  (CFD trace shifted to start at the first model's inlet pressure).
%   'lumped_results' : lumped-model results struct to draw when 'results_2D' is
%                  given directly (otherwise it is read from the results folder).
%   'cfd_coolant' : optional cell array (same length/order as cfd) of
%                  process_wall_CFD_data structs for the COOLANT side (the complete
%                  coolant flow, processed with d_h_coolant and @k_EG). The coolant
%                  exports run along the same stations as the air (z = 0 ... -0.300), so
%                  the coolant figures are drawn against the position along the HX, like
%                  the air ones, for the model's coolant column l_index in every air segment.
%   'cfd_coolant_ends' : optional cell array (same length/order) with the coolant INLET and
%                  OUTLET surface values (4th output of load_cfd_exports). Adds the coolant
%                  inlet / outlet temperature and pressure-drop figures.
%   'cool_flip'  : obsolete (the coolant stations run along the HX like the air);
%                  accepted and ignored.
%   'l_index'    : coolant-side segment (column) compared (default 1, see
%                  compare_cfd_vs_2d for why).
%   'nu_ratio'   : true (default) draws Nu_CFD / Nu_developed against the entrance
%                  factors (see below).
%
% Wall temperature: the CFD air-side wall temperature (area average on each station) is
% compared with the model's wall temperature, which in the model is the coolant
% temperature: T_cool_seg(k, l_index), the coolant entering column l of segment k (the
% value used for the temperature ratio and the DNS lookup), plus the mean over the column
% as a dashed line. 2D models only.
%
% Coolant figures (all against the position along the HX; the model value of air segment k is
% drawn at the segment middle, for coolant column l = l_index):
%   inlet temperature  : CFD coolant bulk temperature on the coolant INLET surface vs the model's
%                        coolant temperature entering column l, T_cool_seg(k,l). This is also the
%                        model's wall temperature (= T_h_i for l = 1).
%   outlet temperature : CFD coolant bulk temperature on the coolant OUTLET surface vs the model's
%                        coolant temperature leaving column l (= entering column l+1),
%                        T_cool_seg(k,l+1). 2D models only.
%   bulk T, v, Nu, h, rho, Cf, Re, y+ : CFD complete coolant flow vs the model's column l
%                        (v_cool_seg, h_cool_seg, Re_cool_seg; Nu = h*d_h/k_EG at the coolant
%                        temperature the model evaluates its properties at; Cf = f/4 since the
%                        model's f is Darcy). The CFD Nu uses q_coolant with the wall and bulk
%                        temperatures. y+ is CFD only.
%   dP                 : CFD inlet-surface minus outlet-surface pressure vs dp_cool_seg(k,l).
%
% Nu ratio figure: for every CFD station the DEVELOPED Nusselt number is evaluated at the
% CFD's own local state (Re, Pr from the local bulk temperature, wall/bulk temperature
% ratio): the model correlation without its length term, and/or the DNS value with
% L = Inf, depending on which model variants are present. The ratio Nu_CFD / Nu_developed
% is the entrance enhancement the CFD actually shows. It is drawn against the factors the
% model uses: the local factor 1+(1/3)(d/x)^(2/3), the segment-average factor over each
% model segment, the factor as evaluated in the code (local form at x = k*dx) and the
% original length-averaged Gnielinski factor 1+(d/x)^(2/3). Only meaningful for Re above
% the laminar range, where the model applies the factor.
%
% Also drawn when the CFD structs carry them: Reynolds number (vs the model's
% Re_air_seg) and wall y+ (CFD only).
%
% Cf: the CFD skin friction coefficient (Fanning, from wall shear stress) is compared with the
% model's air-side f_air_seg, which is already Fanning (Konakov 0.25*(...)^-2, or the DNS Cf).
% On the coolant side the model's f is Darcy, so f/4 is used there. For the 2D model only coolant-side
% segment l_index (default 1) is used. The lumped model gives one constant f.
%
% Pressure: with 'cfd_P_operating' the CFD trace is the actual pressure
% (Pbulk + operating pressure). Without it, each CFD trace is shifted to start
% at the inlet static pressure of the FIRST 2D result, P_air_seg(1,l).
%
% Cumulative heat transfer [W per CFD unit cell], from cumulative_heat_transfer.m:
% CFD from the Fluent UDF file cum_heat_transfer.txt (wall heat flux summed from the inlet,
% read by load_cfd_exports as cfd{k}.Q_cum_wall) when it is present, otherwise by an energy
% balance on the air. Models as the running sum of Q_seg_arr in
% the first coolant column (l_index), scaled down to one cell with the channel
% counts. See that function for the scaling.

p = inputParser;
addParameter(p, 'cfd_labels', {});
addParameter(p, 'results_2D', {});
addParameter(p, 'labels_2D', {});
addParameter(p, 'id', []);
addParameter(p, 'N_segments', []);
addParameter(p, 'N_cool_seg', []);
addParameter(p, 'use_DNS', [false true]);
addParameter(p, 'include_1d', true);
addParameter(p, 'include_lumped', true);
addParameter(p, 'lumped_results', []);
addParameter(p, 'cfd_P_operating', []);
addParameter(p, 'cfd_coolant', {});
addParameter(p, 'cfd_coolant_ends', {});
addParameter(p, 'cool_flip', false);   % obsolete, ignored
addParameter(p, 'l_index', 1);
addParameter(p, 'nu_ratio', true);
parse(p, varargin{:});
l = p.Results.l_index;

if isstruct(cfd), cfd = {cfd}; end
n_cfd = numel(cfd);
cfd_labels = p.Results.cfd_labels;
if isempty(cfd_labels)
    cfd_labels = arrayfun(@(k) sprintf('CFD %d', k), 1:n_cfd, 'UniformOutput', false);
elseif numel(cfd_labels) ~= n_cfd
    error('compare_cfd_meshes_vs_2d:labels', 'cfd_labels has %d entries but %d CFD cases were given.', numel(cfd_labels), n_cfd);
end

%% ---- 2D results: supplied, or read from the results folder ---- %%
res2D  = p.Results.results_2D;
lab2D  = p.Results.labels_2D;
dns2D  = [];
is1d   = [];
lumped = [];
if isstruct(res2D), res2D = {res2D}; end
if ~isempty(res2D)
    if isempty(lab2D)
        lab2D = arrayfun(@(k) sprintf('2D model %d', k), 1:numel(res2D), 'UniformOutput', false);
    elseif numel(lab2D) ~= numel(res2D)
        error('compare_cfd_meshes_vs_2d:labels', 'labels_2D has %d entries but %d 2D results were given.', numel(lab2D), numel(res2D));
    end
    if ~isempty(p.Results.lumped_results), lumped = p.Results.lumped_results; end
    for k = 1:numel(res2D)
        dns2D(k) = isfield(res2D{k}, 'use_DNS') && res2D{k}.use_DNS; %#ok<AGROW>
        is1d(k)  = isfield(res2D{k}, 'N_cool_seg') && res2D{k}.N_cool_seg == 0; %#ok<AGROW>
    end
else
    id = p.Results.id;  N_seg = p.Results.N_segments;  N_cool = p.Results.N_cool_seg;
    if isempty(id) || isempty(N_seg) || isempty(N_cool)
        error('compare_cfd_meshes_vs_2d:noSource', ...
            'Give either ''results_2D'' or all of ''id'', ''N_segments'' and ''N_cool_seg'' to read the 2D results from the results folder.');
    end
    idx_file = fullfile('results', 'results_index.mat');
    if ~isfile(idx_file)
        error('compare_cfd_meshes_vs_2d:noIndex', '%s not found.', idx_file);
    end
    loaded = load(idx_file, 'results_index');
    row = find(loaded.results_index.id == id, 1);
    if isempty(row)
        error('compare_cfd_meshes_vs_2d:noId', 'id %d is not in the results index.', id);
    end
    solve_T = loaded.results_index.inputs{row}.solve_T;
    if p.Results.include_1d, Nc_list = [0 N_cool]; else, Nc_list = N_cool; end
    for Nc = Nc_list
        for dns = p.Results.use_DNS(:)'
            f = fullfile('results', sprintf('id-%d', id), disc_results_filename(N_seg, Nc, dns, id, solve_T));
            if isfile(f)
                ld = load(f, 'results');
                res2D{end+1} = ld.results; %#ok<AGROW>
                dns2D(end+1) = dns;        %#ok<AGROW>
                is1d(end+1)  = (Nc == 0);  %#ok<AGROW>
                if dns, kind = 'DNS-based'; else, kind = 'correlations'; end
                if Nc == 0
                    lab2D{end+1} = sprintf('1D model, %s (N=%d)', kind, N_seg); %#ok<AGROW>
                else
                    lab2D{end+1} = sprintf('2D model, %s (N=%dx%d)', kind, N_seg, Nc); %#ok<AGROW>
                end
                fprintf('Loaded model results: %s\n', f);
            else
                warning('compare_cfd_meshes_vs_2d:missing', 'No saved result for N_cool_seg = %d, use_DNS = %d: %s -- skipped.', Nc, dns, f);
            end
        end
    end
    if p.Results.include_lumped
        fL = fullfile('results', sprintf('id-%d', id), lumped_results_filename(id, solve_T));
        if isfile(fL)
            ld = load(fL, 'results');
            lumped = ld.results;
            fprintf('Loaded lumped results: %s\n', fL);
        else
            warning('compare_cfd_meshes_vs_2d:missingLumped', 'No saved lumped result: %s -- skipped.', fL);
        end
    end
    if isempty(res2D)
        error('compare_cfd_meshes_vs_2d:none', ...
            'No saved 2D result found for id %d, N_segments = %d, N_cool_seg = %d (use_DNS = %s).', ...
            id, N_seg, N_cool, mat2str(p.Results.use_DNS));
    end
end
n_2d = numel(res2D);
% Column used from each result: the 1D model has a single (coolant-lumped) column.
col = zeros(1, n_2d);
for k = 1:n_2d
    if is1d(k), col(k) = 1; else, col(k) = l; end
end
for k = 1:n_2d
    if ~is1d(k) && l > size(res2D{k}.T_air_seg, 2)
        error('compare_cfd_meshes_vs_2d:lOutOfRange', 'l_index = %d but 2D result %d only has %d coolant-side segments.', ...
            l, k, size(res2D{k}.T_air_seg, 2));
    end
end

%% ---- Styles ---- %%
% CFD: greys, light (first/coarsest) to dark (last/finest), solid, with markers.
% Models: colour = correlation (blue) or DNS (red); line style = number of coolant-side
% segments (1D model solid, then dashed, dash-dot, dotted for increasing N_cool_seg).
grey    = linspace(0.75, 0.1, max(n_cfd, 2));
mk_cfd  = {'x', '+', 's', '^', 'v', 'o'};
col_corr_model = [0.122 0.341 0.600];    % muted blue  : correlation-based
col_dns_model  = [0.700 0.220 0.200];    % muted red   : DNS-based
ls_list = {'-', '--', '-.', ':'};
Nc_all  = zeros(1, n_2d);
for k = 1:n_2d
    if is1d(k), Nc_all(k) = 0; else, Nc_all(k) = res2D{k}.N_cool_seg; end
end
Nc_levels = unique(Nc_all);               % sorted: 0 (1D) first
sty2D = cell(1, n_2d);
for k = 1:n_2d
    if dns2D(k), c = col_dns_model; m = 's'; else, c = col_corr_model; m = 'o'; end
    ls = ls_list{mod(find(Nc_levels == Nc_all(k), 1) - 1, numel(ls_list)) + 1};
    sty2D{k} = struct('color', c, 'marker', m, 'ls', ls);
end
lw = 1.2; ms = 4;

P_inlet = res2D{1}.P_air_seg(1, col(1));
P_op    = p.Results.cfd_P_operating;
use_P_op = ~isempty(P_op);
N_fin_ref = res2D{1}.N_fin_air;  N_pass_ref = res2D{1}.N_air_pass;   % channel counts (for the lumped model)

%% ---- Figures ---- %%
fig_T  = new_fig('Air Bulk Temperature: CFD meshes vs 2D model', 'Position along HX [mm]', 'Air bulk temperature [K]', 'Air bulk temperature: CFD meshes vs 2D model');
fig_P  = new_fig('Static Pressure: CFD meshes vs 2D model',      'Position along HX [mm]', 'Static pressure [Pa]',        'Static pressure: CFD meshes vs 2D model');
fig_h  = new_fig('Heat Transfer Coefficient: CFD meshes vs 2D model', 'Position along HX [mm]', 'h [W/m^2K]',           'Heat transfer coefficient: CFD meshes vs 2D model');
fig_v  = new_fig('Velocity: CFD meshes vs 2D model',             'Position along HX [mm]', 'Velocity [m/s]',             'Channel velocity: CFD meshes vs 2D model');
fig_m  = new_fig('Mass Flow Rate: CFD meshes vs 2D model',       'Position along HX [mm]', 'Mass flow rate [kg/s]',      'Mass flow rate: CFD meshes vs 2D model');
has_rho = any(cellfun(@(r) isfield(r, 'rho'), cfd));
if has_rho
    fig_rho = new_fig('Density: CFD meshes vs models', 'Position along HX [mm]', 'Density [kg/m^3]', 'Air density: CFD meshes vs models');
end
has_cf = any(cellfun(@(r) isfield(r, 'Cf'), cfd));
has_Re = any(cellfun(@(r) isfield(r, 'Re'), cfd));
has_yp = any(cellfun(@(r) isfield(r, 'yplus'), cfd));
if has_Re
    fig_Re = new_fig('Reynolds Number: CFD meshes vs models', 'Position along HX [mm]', 'Re [-]', 'Air Reynolds number: CFD meshes vs models');
end
if has_yp
    fig_yp = new_fig('y+: CFD meshes', 'Position along HX [mm]', 'y+ [-]', 'Air-side wall y+ (CFD)');
end
if has_cf
    fig_cf = new_fig('Skin Friction: CFD meshes vs models', 'Position along HX [mm]', 'C_f (Fanning) [-]', 'Air skin friction coefficient: CFD meshes vs models');
end
fig_Q  = new_fig('Cumulative Heat Transfer: CFD meshes vs 2D model', 'Position along HX [mm]', 'Cumulative heat transfer per unit cell [W]', 'Cumulative heat transfer: CFD meshes vs models');
fig_Nu = new_fig('Nusselt Number: CFD meshes vs 2D model',       'Position along HX [mm]', 'Nu [-]',                     'Nusselt number: CFD meshes vs 2D model');
has_2d = any(~is1d);
if has_2d
    fig_Tw = new_fig('Air-side Wall Temperature: CFD vs model', 'Position along HX [mm]', 'Wall temperature [K]', ...
        sprintf('Air-side wall temperature: CFD vs model (= coolant temperature, column %d)', l));
end

comparison.l_index = l;
comparison.cfd_labels = cfd_labels;
comparison.labels_2D = lab2D;

% CFD cases
for k = 1:n_cfd
    r  = cfd{k};
    c  = grey(k)*[1 1 1];
    m  = mk_cfd{mod(k-1, numel(mk_cfd)) + 1};
    x_all = -r.z * 1000;
    valid = ~r.is_entrance;
    x_v   = -r.z(valid) * 1000;
    args  = {'-', 'Color', c, 'LineWidth', lw, 'Marker', m, 'MarkerSize', ms+1, 'DisplayName', cfd_labels{k}};
    plot(fig_T.ax,  x_all, r.Tbulk, args{:});
    if use_P_op, P_cfd = r.Pbulk + P_op; else, P_cfd = P_inlet + r.dP_bulk; end
    plot(fig_P.ax,  x_all, P_cfd, args{:});
    if isfield(r, 'Q_cum_wall')      % Fluent UDF: wall heat flux summed from the inlet (cum_heat_transfer.txt)
        x_Q = -r.Q_cum_wall.z * 1000;  Q_cfd = r.Q_cum_wall.Q;
    else                             % fallback: energy balance on the air
        x_Q = x_all;  Q_cfd = cumulative_heat_transfer(r, -r.z);
    end
    plot(fig_Q.ax,  x_Q, Q_cfd, args{:});
    plot(fig_h.ax,  x_v,   r.h(valid),  args{:});
    plot(fig_v.ax,  x_all, r.Vmag, args{:});
    plot(fig_m.ax,  x_all, r.mdot, args{:});
    plot(fig_Nu.ax, x_v,   r.Nu(valid), args{:});
    if has_2d, plot(fig_Tw.ax, x_all, r.Twall, args{:}); end
    if isfield(r, 'rho'), plot(fig_rho.ax, x_all, r.rho, args{:}); end
    if isfield(r, 'Cf'),  plot(fig_cf.ax,  x_all, r.Cf,  args{:}); end
    if isfield(r, 'Re'),  plot(fig_Re.ax,  x_all, r.Re,  args{:}); end
    if isfield(r, 'yplus'), plot(fig_yp.ax, x_all, r.yplus, args{:}); end
    comparison.cfd(k).label = cfd_labels{k};
    comparison.cfd(k).x_all = x_all;  comparison.cfd(k).T = r.Tbulk;
    comparison.cfd(k).P = P_cfd;  comparison.cfd(k).Q_cum = Q_cfd;  comparison.cfd(k).V = r.Vmag;  comparison.cfd(k).mdot = r.mdot;
    comparison.cfd(k).x_valid = x_v;  comparison.cfd(k).h = r.h(valid);  comparison.cfd(k).Nu = r.Nu(valid);
end

% 2D model variants
for k = 1:n_2d
    r = res2D{k};
    x_b = linspace(0, r.L_solution, r.N_segments+1) * 1000;
    x_m = 0.5*(x_b(1:end-1) + x_b(2:end));
    mdot_model = r.m_dot_seg./r.N_fin_air*2./r.N_air_pass;
    s = sty2D{k};
    args_node = {s.ls, 'Color', s.color, 'LineWidth', lw, 'Marker', s.marker, 'MarkerFaceColor', s.color, 'MarkerSize', ms, 'DisplayName', lab2D{k}};
    plot(fig_T.ax,  x_b, r.T_air_seg(:, col(k))',       args_node{:});
    plot(fig_P.ax,  x_b, r.P_air_seg(:, col(k))',       args_node{:});
    plot(fig_h.ax,  x_m, r.h_air_seg(:, col(k))',       args_node{:});
    plot(fig_v.ax,  x_b, r.v_channel_seg(:, col(k))',   args_node{:});
    plot(fig_m.ax,  x_b, mdot_model,               args_node{:});
    plot(fig_Nu.ax, x_m, r.Nu_air_seg(:, col(k))',      args_node{:});
    if ~is1d(k) && isfield(r, 'T_cool_seg')
        % Model wall temperature = coolant temperature: the coolant entering column l of
        % each air segment (what the model uses for the temperature ratio / DNS lookup).
        Tw_model = r.T_cool_seg(:, col(k))';
        args_w = args_node;  args_w{end} = sprintf('%s, coolant inlet T of column %d', lab2D{k}, col(k));
        plot(fig_Tw.ax, x_m, Tw_model, args_w{:});
        comparison.model2D(k).T_wall = Tw_model;
        if isfield(r, 'T_mean_h_arr')
            plot(fig_Tw.ax, x_m, r.T_mean_h_arr(:, col(k))', s.ls, 'Color', 0.55*s.color + 0.45, 'LineWidth', 1, ...
                'DisplayName', sprintf('%s, column %d mean', lab2D{k}, col(k)));
        end
    end
    if has_rho
        rho_model = arrayfun(@(Pk, Tk) rho_air(Pk, Tk), r.P_air_seg(:, col(k)), r.T_air_seg(:, col(k)))';
        plot(fig_rho.ax, x_b, rho_model, args_node{:});
    end
    Q_model = cumulative_heat_transfer(r, x_b'/1000, 'column', col(k))';
    plot(fig_Q.ax, x_b, Q_model, args_node{:});
    comparison.model2D(k).Q_cum = Q_model;
    if has_Re
        plot(fig_Re.ax, x_b, r.Re_air_seg(:, col(k))', args_node{:});
    end
    if has_cf
        plot(fig_cf.ax, x_b, r.f_air_seg(:, col(k))', args_node{:});
        comparison.model2D(k).Cf = r.f_air_seg(:, col(k))';
    end
    comparison.model2D(k).label = lab2D{k};
    comparison.model2D(k).x_bounds = x_b;  comparison.model2D(k).x_mid = x_m;
    comparison.model2D(k).T = r.T_air_seg(:, col(k))';  comparison.model2D(k).P = r.P_air_seg(:, col(k))';
    comparison.model2D(k).h = r.h_air_seg(:, col(k))';  comparison.model2D(k).Nu = r.Nu_air_seg(:, col(k))';
    comparison.model2D(k).v = r.v_channel_seg(:, col(k))';  comparison.model2D(k).mdot = mdot_model;
end

% Lumped model: straight lines inlet -> outlet (no position resolution)
if ~isempty(lumped)
    xL = [0 lumped.L_solution]*1000;
    argsL = {'-', 'Color', [0 0 0], 'LineWidth', 1.5, 'Marker', 's', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 6, ...
             'DisplayName', sprintf('Lumped model (L=%.3fm)', lumped.L_solution)};
    plot(fig_T.ax,  xL, [lumped.T3 lumped.T4], argsL{:});
    plot(fig_P.ax,  xL, [lumped.P3 lumped.P4], argsL{:});
    plot(fig_h.ax,  xL, [lumped.h_air lumped.h_air], argsL{:});
    plot(fig_v.ax,  xL, [lumped.v_channel_air lumped.v_channel_air], argsL{:});
    plot(fig_Nu.ax, xL, [lumped.Nu_air lumped.Nu_air], argsL{:});
    if has_rho
        plot(fig_rho.ax, xL, [rho_air(lumped.P3, lumped.T3) rho_air(lumped.P4, lumped.T4)], argsL{:});
    end
    if isfield(lumped, 'Q_pred_solution')
        plot(fig_Q.ax, xL, cumulative_heat_transfer(lumped, xL/1000, 'N_fin_air', N_fin_ref, 'N_air_pass', N_pass_ref), argsL{:});
    end
    if has_Re && isfield(lumped, 'Re_air')
        plot(fig_Re.ax, xL, [lumped.Re_air lumped.Re_air], argsL{:});
    end
    if has_cf && isfield(lumped, 'f_air')
        plot(fig_cf.ax, xL, [lumped.f_air lumped.f_air], argsL{:});
    end
    comparison.lumped = lumped;
end

for f = [fig_T fig_P fig_h fig_v fig_m fig_Nu fig_Q]
    legend_unique(f.ax);
end
if has_2d,  legend_unique(fig_Tw.ax); end
if has_rho, legend_unique(fig_rho.ax); end
if has_cf,  legend_unique(fig_cf.ax); end
if has_Re,  legend_unique(fig_Re.ax); end
if has_yp,  legend_unique(fig_yp.ax); end

comparison.results_2D = res2D;  comparison.lumped_results = lumped;

%% ---- Nu(CFD) / developed Nu at the CFD state, against the entrance factors ---- %%
if p.Results.nu_ratio
    comparison.nu_ratio = plot_nu_ratio(cfd, cfd_labels, grey, mk_cfd, res2D, dns2D, P_op, P_inlet, lw, ms);
end

%% ---- Coolant side (against the position along the HX) ---- %%
cfd_cool = p.Results.cfd_coolant;
if isstruct(cfd_cool), cfd_cool = {cfd_cool}; end
cfd_ends = p.Results.cfd_coolant_ends;
if isstruct(cfd_ends), cfd_ends = {cfd_ends}; end
if ~isempty(cfd_cool) || ~isempty(cfd_ends)
    if ~isempty(cfd_cool) && numel(cfd_cool) ~= n_cfd
        error('compare_cfd_meshes_vs_2d:coolCount', 'cfd_coolant has %d entries but %d air-side CFD cases were given.', numel(cfd_cool), n_cfd);
    end
    if ~isempty(cfd_ends) && numel(cfd_ends) ~= n_cfd
        error('compare_cfd_meshes_vs_2d:coolCount', 'cfd_coolant_ends has %d entries but %d air-side CFD cases were given.', numel(cfd_ends), n_cfd);
    end
    comparison.coolant = plot_coolant_hx(cfd_cool, cfd_ends, cfd_labels, res2D, lab2D, is1d, sty2D, col, grey, mk_cfd, lw, ms);
end

end

%% ---- Local helper: one figure with a single labelled axes ---- %%
function f = new_fig(name, xl, yl, ttl)
f.fig = figure('Name', name, 'NumberTitle', 'off');
f.ax  = axes(f.fig); hold(f.ax, 'on'); grid(f.ax, 'on');
xlabel(f.ax, xl); ylabel(f.ax, yl); title(f.ax, ttl);
end

%% ---- Local helper: coolant-side figures against the position along the HX ---- %%
% The coolant exports run along the same stations as the air (z = 0 ... -0.300), so every
% figure shares the air figures' x axis. Model: air segment k (drawn at the segment middle),
% coolant column l = col(k) (1 for the 1D model, which has a single column).
%   cfd_c    : complete-flow coolant structs (process_wall_CFD_data), {} if not given
%   cfd_ends : coolant inlet / outlet surface values per mesh (load_cfd_exports 4th output), {} if not given
function out = plot_coolant_hx(cfd_c, cfd_ends, cfd_labels, res2D, lab2D, is1d, sty2D, col, grey, mk_cfd, lw, ms)
out = struct();
have_c    = ~isempty(cfd_c);
have_ends = ~isempty(cfd_ends) && any(~cellfun(@isempty, cfd_ends));
n_cfd     = max(numel(cfd_c), numel(cfd_ends));
l_title   = 1;  if any(~is1d), l_title = max(col(~is1d)); end
xl = 'Position along HX [mm]';

F = struct();
if have_c
    F.T  = new_fig('Coolant Bulk Temperature: CFD vs model', xl, 'Coolant bulk temperature [K]', ...
        'Coolant bulk temperature (complete flow): CFD vs model column mean');
    F.v  = new_fig('Coolant Velocity: CFD vs model', xl, 'Velocity [m/s]', 'Coolant velocity: CFD vs model');
    F.h  = new_fig('Coolant Heat Transfer Coefficient: CFD vs model', xl, 'h [W/m^2K]', 'Coolant heat transfer coefficient: CFD (q_{coolant}) vs model');
    F.Nu = new_fig('Coolant Nusselt Number: CFD vs model', xl, 'Nu [-]', 'Coolant Nusselt number: CFD (q_{coolant}) vs model');
    if any(cellfun(@(r) isfield(r, 'rho'), cfd_c))
        F.rho = new_fig('Coolant Density: CFD vs model', xl, 'Density [kg/m^3]', 'Coolant density: CFD vs model');
    end
    if any(cellfun(@(r) isfield(r, 'Cf'), cfd_c))
        F.Cf = new_fig('Coolant Skin Friction: CFD vs model', xl, 'C_f = f/4 [-]', 'Coolant skin friction coefficient: CFD vs model');
    end
    if any(cellfun(@(r) isfield(r, 'Re'), cfd_c))
        F.Re = new_fig('Coolant Reynolds Number: CFD vs model', xl, 'Re [-]', 'Coolant Reynolds number: CFD vs model');
    end
    if any(cellfun(@(r) isfield(r, 'yplus'), cfd_c))
        F.yplus = new_fig('Coolant y+: CFD', xl, 'y+ [-]', 'Coolant wall y+ (CFD)');
    end
end
if have_ends
    F.Tin  = new_fig('Coolant Inlet Temperature: CFD vs model wall temperature', xl, 'Coolant inlet temperature [K]', ...
        sprintf('Coolant inlet temperature: CFD inlet surface vs model coolant entering column %d (= model wall temperature)', l_title));
    F.Tout = new_fig('Coolant Outlet Temperature: CFD vs model', xl, 'Coolant outlet temperature [K]', ...
        sprintf('Coolant outlet temperature: CFD outlet surface vs model coolant entering column %d (= leaving column %d)', l_title+1, l_title));
    if any(cellfun(@(e) isstruct(e) && isfield(e, 'dP'), cfd_ends))
        F.dP = new_fig('Coolant Pressure Drop: CFD vs model', xl, 'Pressure drop [Pa]', ...
            sprintf('Coolant pressure drop over one coolant column: CFD (inlet - outlet surface) vs model column %d', l_title));
    end
end

% ---- CFD ----
for k = 1:n_cfd
    c = grey(k)*[1 1 1];  m = mk_cfd{mod(k-1, numel(mk_cfd)) + 1};
    args = {'-', 'Color', c, 'LineWidth', lw, 'Marker', m, 'MarkerSize', ms+1, 'DisplayName', cfd_labels{k}};
    if have_c
        r = cfd_c{k};  x_all = -r.z*1000;  valid = ~r.is_entrance;  x_v = x_all(valid);
        plot(F.T.ax,  x_all, r.Tbulk, args{:});
        plot(F.v.ax,  x_all, r.Vmag,  args{:});
        plot(F.h.ax,  x_v, r.h(valid),  args{:});
        plot(F.Nu.ax, x_v, r.Nu(valid), args{:});
        if isfield(F, 'rho')   && isfield(r, 'rho'),   plot(F.rho.ax,   x_all, r.rho,   args{:}); end
        if isfield(F, 'Cf')    && isfield(r, 'Cf'),    plot(F.Cf.ax,    x_all, r.Cf,    args{:}); end
        if isfield(F, 'Re')    && isfield(r, 'Re'),    plot(F.Re.ax,    x_all, r.Re,    args{:}); end
        if isfield(F, 'yplus') && isfield(r, 'yplus'), plot(F.yplus.ax, x_all, r.yplus, args{:}); end
    end
    if have_ends && ~isempty(cfd_ends{k})
        e = cfd_ends{k};
        plot(F.Tin.ax,  -e.z_in*1000,  e.Tin,  args{:});
        plot(F.Tout.ax, -e.z_out*1000, e.Tout, args{:});
        if isfield(F, 'dP') && isfield(e, 'dP'), plot(F.dP.ax, -e.z_dP*1000, e.dP, args{:}); end
    end
end

% ---- Model: air segment k, coolant column l ----
need = {'T_cool_seg', 'v_cool_seg', 'Re_cool_seg', 'h_cool_seg', 'd_h_coolant', 'dp_cool_seg'};
for k = 1:numel(res2D)
    r = res2D{k};
    if ~all(isfield(r, need)) || (is1d(k) && ~isfield(r, 'T_mean_h_arr'))
        warning('compare_cfd_meshes_vs_2d:noCoolantFields', ...
            '%s lacks saved coolant fields (%s) -- rerun it; its coolant curves are skipped.', lab2D{k}, strjoin(need(~isfield(r, need)), ', '));
        continue
    end
    lk = col(k);
    xb = linspace(0, r.L_solution, r.N_segments+1)*1000;
    xm = 0.5*(xb(1:end-1) + xb(2:end));
    T_in = r.T_cool_seg(:, lk)';                       % coolant entering column l of each air segment
    if is1d(k)
        T_ref  = r.T_mean_h_arr(:, 1)';                % 1D: properties at the iterated mean coolant temperature
        T_out  = nan(size(T_in));                      % no column resolution in the 1D model
        T_mean = T_ref;
    else
        T_ref  = T_in;                                 % 2D: properties at the coolant entering the column
        T_out  = r.T_cool_seg(:, lk+1)';               % coolant leaving column l (= entering column l+1)
        T_mean = 0.5*(T_in + T_out);
    end
    dh  = r.d_h_coolant;
    Re  = r.Re_cool_seg(:, lk)';
    f_D = (1 ./ (-1.8*log10((0.0015/dh)^1.11 + 6.9./Re))).^2;
    v   = r.v_cool_seg(:, lk)';
    h   = r.h_cool_seg(:, lk)';
    Nu  = h*dh ./ arrayfun(@k_EG, T_ref);
    dpc = r.dp_cool_seg(:, lk)';

    st = sty2D{k};
    a_n = {st.ls, 'Color', st.color, 'LineWidth', lw, 'Marker', st.marker, 'MarkerFaceColor', st.color, 'MarkerSize', ms, 'DisplayName', lab2D{k}};
    if isfield(F, 'Tin')
        a_i = a_n;  a_i{end} = sprintf('%s: entering column %d', lab2D{k}, lk);
        plot(F.Tin.ax, xm, T_in, a_i{:});
    end
    if isfield(F, 'Tout') && ~is1d(k)
        a_o = a_n;  a_o{end} = sprintf('%s: leaving column %d (= entering column %d)', lab2D{k}, lk, lk+1);
        plot(F.Tout.ax, xm, T_out, a_o{:});
    end
    if isfield(F, 'dP'), plot(F.dP.ax, xm, dpc, a_n{:}); end
    if have_c
        plot(F.T.ax,  xm, T_mean, a_n{:});
        plot(F.v.ax,  xm, v,      a_n{:});
        plot(F.h.ax,  xm, h,      a_n{:});
        plot(F.Nu.ax, xm, Nu,     a_n{:});
        if isfield(F, 'rho'), plot(F.rho.ax, xm, arrayfun(@rho_EG, T_ref), a_n{:}); end
        if isfield(F, 'Cf'),  plot(F.Cf.ax,  xm, f_D/4, a_n{:}); end
        if isfield(F, 'Re'),  plot(F.Re.ax,  xm, Re,     a_n{:}); end
    end
    out.model(k).label = lab2D{k};  out.model(k).x_mid = xm;  out.model(k).column = lk;
    out.model(k).T_in = T_in;  out.model(k).T_out = T_out;  out.model(k).T_mean = T_mean;
    out.model(k).v = v;  out.model(k).h = h;  out.model(k).Nu = Nu;  out.model(k).Re = Re;  out.model(k).dP = dpc;
end

for fn = fieldnames(F)'
    legend_unique(F.(fn{1}).ax);
end
end

%% ---- Local helper: Nu(CFD) / developed Nu at the CFD state, against the entrance factors ---- %%
% The developed Nu is evaluated at each CFD station's own state: Re (CFD), Pr from the local
% bulk temperature, temperature ratio Twall/Tbulk (the model uses T_cool/T_air). For the
% correlation variant that is the model correlation without the length term; for the DNS
% variant it is the thermoturb lookup with L = Inf. The ratio is the entrance enhancement
% that the CFD actually shows.
function out = plot_nu_ratio(cfd, cfd_labels, grey, mk_cfd, res2D, dns2D, P_op, P_inlet, lw, ms)
out = struct();
r1 = res2D{1};
if ~all(isfield(r1, {'d_h_air', 'L_solution', 'N_segments'}))
    warning('compare_cfd_meshes_vs_2d:noDh', ...
        'Nu ratio figure skipped: the first model result has no d_h_air / L_solution / N_segments.');
    return
end
d_h = r1.d_h_air;  L = r1.L_solution;  N = r1.N_segments;

variants = {};
if any(~dns2D), variants{end+1} = 'corr'; end
if any(dns2D),  variants{end+1} = 'dns';  end
if isempty(variants), return, end
vtitle = struct('corr', 'Developed Nu: model correlation without its length term', ...
                'dns',  'Developed Nu: DNS data (L = Inf)');

% Entrance factors of the model, per segment
dx = L/N;  kk = (1:N)';
x1 = (kk-1)*dx;  x2 = kk*dx;  x_mid = 0.5*(x1 + x2);
F_seg  = 1 + d_h^(2/3) * (x2.^(1/3) - x1.^(1/3)) ./ (x2 - x1);   % local factor averaged over the segment
F_code = 1 + (1/3) * (d_h ./ x2).^(2/3);                          % local form at x = k*dx, as in the code

% CFD state at every station (entrance station left out)
n_cfd = numel(cfd);
data = struct('x', {}, 'Re', {}, 'Pr', {}, 'Tb', {}, 'Tw', {}, 'Nu', {});
for k = 1:n_cfd
    r = cfd{k};
    valid = ~r.is_entrance;
    Tb = r.Tbulk(valid);  Tw = r.Twall(valid);
    if isfield(r, 'Re')
        Re = r.Re(valid);
    else
        if isempty(P_op), P = P_inlet + r.dP_bulk(valid); else, P = r.Pbulk(valid) + P_op; end
        Re = arrayfun(@(Pk, Tk, Vk) rho_air(Pk, Tk)*Vk*d_h/mu_air(Tk), P, Tb, r.Vmag(valid));
    end
    data(k).x  = -r.z(valid);
    data(k).Re = Re;
    data(k).Pr = arrayfun(@(Tk) mu_air(Tk)*cp_air(Tk)/k_air(Tk), Tb);
    data(k).Tb = Tb;  data(k).Tw = Tw;
    data(k).Nu = r.Nu(valid);
    fprintf('Nu ratio, %s: Re %.0f to %.0f, Tw/Tb %.2f to %.2f, %.0f%% of stations below Re 3000 (no length factor there)\n', ...
        cfd_labels{k}, min(Re), max(Re), min(Tw./Tb), max(Tw./Tb), 100*mean(Re < 3000));
end

f = figure('Name', 'Nu ratio: CFD vs developed Nu at the CFD state', 'NumberTitle', 'off');
tl = tiledlayout(f, numel(variants), 1, 'TileSpacing', 'compact');
title(tl, 'CFD Nu / developed Nu at the CFD state, against the model entrance factors');
xmax = max(arrayfun(@(d) max(d.x), data));
xmin = min(arrayfun(@(d) min(d.x), data));
x_line = linspace(xmin, xmax, 400);

for iv = 1:numel(variants)
    v = variants{iv};
    a = nexttile(tl); hold(a, 'on'); grid(a, 'on');
    xlabel(a, 'Distance from the inlet [mm]'); ylabel(a, 'Nu_{CFD} / Nu_{developed} [-]');
    title(a, vtitle.(v));
    out.(v) = struct('label', {}, 'x', {}, 'ratio', {});
    for k = 1:n_cfd
        d = data(k);
        if strcmp(v, 'corr')
            Nu_fd = arrayfun(@nu_developed_corr, d.Re, d.Pr, d.Tw./d.Tb);
        else
            Nu_fd = nan(size(d.Re));
            for i = 1:numel(Nu_fd)
                Nu_fd(i) = nu_developed_dns(d.Re(i), d.Pr(i), d.Tb(i), d.Tw(i), d_h);
            end
            if any(isnan(Nu_fd))
                warning('compare_cfd_meshes_vs_2d:dnsLookup', ...
                    '%s: %d of %d DNS lookups failed or fell outside the table -- those points are left out.', ...
                    cfd_labels{k}, nnz(isnan(Nu_fd)), numel(Nu_fd));
            end
        end
        ratio = d.Nu ./ Nu_fd;
        c = grey(k)*[1 1 1];  m = mk_cfd{mod(k-1, numel(mk_cfd)) + 1};
        plot(a, d.x*1000, ratio, '-', 'Color', c, 'LineWidth', lw, 'Marker', m, 'MarkerSize', ms+1, 'DisplayName', cfd_labels{k});
        out.(v)(k).label = cfd_labels{k};  out.(v)(k).x = d.x;  out.(v)(k).ratio = ratio;
    end
    plot(a, x_line*1000, 1 + (1/3)*(d_h./x_line).^(2/3), 'k-', 'LineWidth', 1.5, 'DisplayName', 'local factor 1+(1/3)(d/x)^{2/3}');
    plot(a, x_line*1000, 1 + (d_h./x_line).^(2/3), 'k:', 'LineWidth', 1.2, 'DisplayName', 'length-averaged factor 1+(d/x)^{2/3}');
    plot(a, x_mid*1000, F_seg, 'o', 'Color', [0.122 0.341 0.600], 'MarkerFaceColor', [0.122 0.341 0.600], 'MarkerSize', ms, ...
        'DisplayName', 'segment-average factor');
    plot(a, x_mid*1000, F_code, '^', 'Color', [0.700 0.220 0.200], 'MarkerFaceColor', [0.700 0.220 0.200], 'MarkerSize', ms, ...
        'DisplayName', 'local factor at x = k dx (as in the code)');
    xlim(a, [0, xmax*1000]);
    legend_unique(a);
end
out.d_h = d_h;  out.x_mid = x_mid;  out.F_seg = F_seg;  out.F_code = F_code;
end

%% ---- Local helper: developed Nu from the model correlation (calculate_air_side, no length term) ---- %%
function Nu = nu_developed_corr(Re, Pr, K)
f = 0.25*((1.8*log10(Re)) - 1.5)^(-2);        % Konakov et al.
if Re < 2300
    Nu = 3.66;
else
    Nu_t = (f/2)*(Re - 1000)*Pr*K / (1 + 12.7*sqrt(f/2)*(Pr^(2/3) - 1));
    if Re < 3000
        w  = (Re - 2300)/(3000 - 2300);
        Nu = (1-w)*3.66 + w*Nu_t;
    else
        Nu = Nu_t;
    end
end
end

%% ---- Local helper: developed Nu from the DNS data (L = Inf, i.e. no length factor) ---- %%
function Nu = nu_developed_dns(Re, Pr, Tb, Tw, d_h)
try
    [~, Nu, ~] = thermoturb_cached(Re, Pr, Tb, Tw, d_h, Inf);
catch
    Nu = NaN;
end
end

%% ---- Local helper: legend with one entry per series name ---- %%
function legend_unique(ax)
ch = flipud(ax.Children(:));                       % plotting order
ch = ch(arrayfun(@(o) isprop(o, 'DisplayName') && ~isempty(o.DisplayName), ch));
if isempty(ch), return, end
[~, ia] = unique({ch.DisplayName}, 'stable');
legend(ax, ch(ia), 'Location', 'best');
end