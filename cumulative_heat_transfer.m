function [Q_cum, z_out] = cumulative_heat_transfer(src, z, varargin)
%CUMULATIVE_HEAT_TRANSFER  Cumulative heat transfer rate [W] from the air-channel
%inlet to position z, for ONE CFD-sized unit cell. Works on a CFD struct
%(from process_wall_CFD_data) or on a saved model result, so both are
%computed the same way and can be plotted against each other.
%
%   Q = cumulative_heat_transfer(src, z)
%   [Q, z_nodes] = cumulative_heat_transfer(src, [])      % at the source's own stations
%
%   src : CFD struct (has Tbulk, mdot, z)      -> energy balance on the air:
%             Q(z) = mdot_inlet * integral of cp_air dT_bulk from the inlet to z
%             (trapezoid rule on the exported stations).
%         disc. model result (has Q_seg_arr)   -> cumulative sum of the segment
%             heat transfer Q_seg_arr(:, column), scaled to one cell (below).
%         lumped result (has Q_pred_solution, L_solution) -> straight line from
%             0 to the total Q over the HX length, scaled the same way.
%   z   : distance from the air inlet along the channel [m] (vector, >= 0).
%         Linear interpolation between stations; [] = the source's own stations.
%
% Name-value options (model / lumped only):
%   'column'       : coolant-side column used (default 1). The 1D model has one
%                    column only.
%   'cell_channels': air channels in one CFD cell (default 2, the same factor the
%                    mass-flow comparison uses: m_dot/N_fin_air*2/N_air_pass).
%   'N_fin_air', 'N_air_pass' : override the channel counts (a lumped result does
%                    not store them; pass them from a discretised result).
%
% Scaling of the model to one cell. Q_seg_arr(k,l) is the heat transferred in
% segment k by the air in coolant column l, i.e. by 1/N_cool of the
% N_fin_air*N_air_pass air channels (the 1D model: all of them). So per air
% channel the heat is  Q_seg_arr(k,l) / (N_fin_air*N_air_pass / N_cool),  and
% one cell holds 'cell_channels' channels. The coolant-side counts
% (N_fin_coolant, N_coolant_pass) do not enter as long as every air channel
% exchanges the same heat with its neighbouring coolant channels.

p = inputParser;
addParameter(p, 'column', 1);
addParameter(p, 'cell_channels', 2);
addParameter(p, 'N_fin_air', []);
addParameter(p, 'N_air_pass', []);
parse(p, varargin{:});

if isfield(src, 'Tbulk') && isfield(src, 'mdot')
    % ---- CFD: energy balance on the air ----
    s_nodes = -src.z(:);                          % distance from the inlet, increasing
    cp_b    = arrayfun(@cp_air, src.Tbulk(:));
    Q_nodes = src.mdot(1) * [0; cumsum(0.5*(cp_b(1:end-1) + cp_b(2:end)) .* diff(src.Tbulk(:)))];

elseif isfield(src, 'Q_seg_arr')
    % ---- discretised model ----
    N_seg   = size(src.Q_seg_arr, 1);
    N_cool  = size(src.Q_seg_arr, 2);
    col     = min(p.Results.column, N_cool);
    [nf, np] = channel_counts(src, p);
    per_cell = p.Results.cell_channels / (nf*np/N_cool);
    s_nodes = linspace(0, src.L_solution, N_seg+1)';
    Q_nodes = [0; cumsum(src.Q_seg_arr(:, col))] * per_cell;

elseif isfield(src, 'Q_pred_solution') && isfield(src, 'L_solution')
    % ---- lumped model ----
    [nf, np] = channel_counts(src, p);
    s_nodes = [0; src.L_solution];
    Q_nodes = [0; src.Q_pred_solution] * p.Results.cell_channels / (nf*np);

else
    error('cumulative_heat_transfer:badSource', ...
        'src must be a CFD struct (Tbulk, mdot, z) or a model result (Q_seg_arr, or Q_pred_solution + L_solution).');
end

if isempty(z)
    z_out = s_nodes;  Q_cum = Q_nodes;
else
    z_out = z;
    Q_cum = interp1(s_nodes, Q_nodes, z, 'linear', 'extrap');
end
end

function [nf, np] = channel_counts(src, p)
nf = p.Results.N_fin_air;   np = p.Results.N_air_pass;
if isempty(nf) && isfield(src, 'N_fin_air'),  nf = src.N_fin_air;  end
if isempty(np) && isfield(src, 'N_air_pass'), np = src.N_air_pass; end
if isempty(nf) || isempty(np)
    error('cumulative_heat_transfer:noChannels', ...
        'N_fin_air / N_air_pass are not in this result -- pass them as name-value options.');
end
end