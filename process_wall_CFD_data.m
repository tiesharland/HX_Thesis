function results = process_wall_CFD_data(q_file, Twall_file, Tbulk_file, Pbulk_file, Vmag_file, mdot_file, d_h_air)
% process_wall_CFD_data  Turn Fluent Surface Integral Report exports
% (q_air, Twall_air, Tbulk_air, Pbulk, Vmag_air, mass flow rate vs a
% shared z-station grid) into a station-by-station h(z) / Nu(z) table
% plus overall averages, an inlet-referenced static pressure trace, bulk
% velocity, and mass flow rate, comparable to the discretised model's
% h_air_seg / Nu_air_seg / P_air_seg / v_channel_seg / m_dot_seg.
%
% Usage:
%   results = process_wall_CFD_data('4by4_150ch_wall_q_air', ...
%                                    '4by4_150ch_wall_Twall_air', ...
%                                    '4by4_150ch_wall_Tbulk_air', ...
%                                    '4by4_150ch_wall_Pbulk', ...
%                                    '4by4_150ch_wall_Vmag_air', ...
%                                    '4by4_150ch_wall_mdot_air', d_h_air);
%   (mdot_file name assumed to follow the same convention as the others -
%   adjust if your actual export uses a different name.)
%
% d_h_air: hydraulic diameter [m] - pass the same value your discretised
% model uses, so Nu(z) is directly comparable.
%
% Requires k_air(T) to be on the MATLAB path (same one used elsewhere in
% the project) for the Nusselt number calculation.
%
% Pressure: Fluent's Pbulk report is a GAUGE pressure relative to a fixed
% operating pressure (ambient at the channel outlet) - not directly
% comparable to the model's absolute static pressure. Since it's just a
% constant offset, shifting the CFD trace so its z=0 (inlet) station
% reads zero removes that offset automatically, whatever its actual
% value - no need to know the operating pressure itself. results.dP_bulk
% is that shifted, inlet-referenced trace (always <= 0, since pressure
% drops downstream); adding it to the model's own inlet static pressure
% (e.g. P_air_seg(1,:)) gives an absolute CFD pressure trace directly
% comparable to the model's P_air_seg.
%
% Mass flow: Fluent reports flux through a +z-facing surface normal, so
% raw values come out negative for flow in -z; results.mdot is the
% absolute value, matching the model's positive m_dot_seg convention.
%
% Assumes all six files report the SAME z-station grid (as your current
% Fluent export does - e.g. z = 0 : -0.004 : -0.300). If a file uses a
% different surface-name prefix per line (e.g. "z-coordinate-547=..." or
% "air-flow-bulk=..."), that's handled automatically; the prefix itself
% is ignored.

[z_q, q]      = read_fluent_report(q_file);
[z_Tw, Twall] = read_fluent_report(Twall_file);
[z_tb, Tbulk] = read_fluent_report(Tbulk_file);
[z_pb, Pbulk] = read_fluent_report(Pbulk_file);
[z_v, Vmag]   = read_fluent_report(Vmag_file);
[z_m, mdot_raw] = read_fluent_report(mdot_file);

if ~isequal(z_q, z_Tw) || ~isequal(z_q, z_tb) || ~isequal(z_q, z_pb) || ...
        ~isequal(z_q, z_v) || ~isequal(z_q, z_m)
    error('process_wall_CFD_data:zMismatch', ...
        ['q, Twall, Tbulk, Pbulk, Vmag and mass-flow files report different ' ...
         'z-stations - check all six exports were built from the same z-list.']);
end
z = z_q;

dT = Twall - Tbulk;
h  = q ./ dT;
Nu = h * d_h_air ./ arrayfun(@k_air, Tbulk);

% Inlet-referenced static pressure: shift so z=0 (inlet, first row since
% read_fluent_report sorts descending) reads zero, cancelling Fluent's
% operating-pressure offset regardless of its value.
dP_bulk = Pbulk - Pbulk(1);

mdot = abs(mdot_raw);   % sign flip only (flow runs in -z); magnitude is the physical quantity

fprintf('%10s %12s %10s %10s %8s %10s %8s %12s %10s %14s\n', ...
    'z [m]', 'q [W/m2]', 'Twall [K]', 'Tbulk [K]', 'dT [K]', 'h [W/m2K]', 'Nu [-]', 'dP [Pa]', 'V [m/s]', 'mdot [kg/s]');
for i = 1:numel(z)
    fprintf('%10.4f %12.2f %10.3f %10.3f %8.3f %10.2f %8.2f %12.3f %10.3f %14.6e\n', ...
        z(i), q(i), Twall(i), Tbulk(i), dT(i), h(i), Nu(i), dP_bulk(i), Vmag(i), mdot(i));
end

% The z=0 (leading-edge) station is a known thermal-entrance singularity:
% h formally diverges there for a step change in wall condition, so its
% CFD value is strongly mesh-size dependent rather than a converged
% physical result. This affects h/Nu ONLY - Tbulk, Pbulk, Vmag and mdot
% are all perfectly well-behaved (non-singular) at z=0 and keep that
% station.
is_entrance = (z == 0);

h_bar_simple  = mean(h(~is_entrance));
Nu_bar_simple = mean(Nu(~is_entrance));

Tw_mean  = mean(Twall(~is_entrance));
T_bulk_in  = Tbulk(is_entrance);       if isempty(T_bulk_in),  T_bulk_in  = Tbulk(1); end
T_bulk_out = Tbulk(end);
dT_in    = Tw_mean - T_bulk_in;
dT_out   = Tw_mean - T_bulk_out;
dT_lm    = (dT_in - dT_out) / log(dT_in/dT_out);
q_bar    = mean(q(~is_entrance));       % or substitute Fluent's reported "Net" value directly
h_bar_lmtd  = q_bar / dT_lm;
Nu_bar_lmtd = h_bar_lmtd * d_h_air / k_air(mean(Tbulk(~is_entrance)));

fprintf('\n--- Summary (z=0 entrance station excluded from h/Nu averages only) ---\n');
fprintf('z=0 entrance-station h         = %.2f W/m2K  (expect mesh-dependent, informational only)\n', h(is_entrance));
fprintf('Simple mean h(z)               = %.2f W/m2K   (min %.2f, max %.2f)\n', ...
    h_bar_simple, min(h(~is_entrance)), max(h(~is_entrance)));
fprintf('Simple mean Nu(z)              = %.2f\n', Nu_bar_simple);
fprintf('Flow-consistent h (q/dT_lm)    = %.2f W/m2K   (dT_lm = %.2f K)\n', h_bar_lmtd, dT_lm);
fprintf('Flow-consistent Nu             = %.2f\n', Nu_bar_lmtd);
fprintf('Total inlet-to-outlet dP (CFD) = %.3f Pa\n', dP_bulk(end));
fprintf('Inlet bulk temperature (z=0)   = %.3f K\n', T_bulk_in);
fprintf('Inlet bulk velocity (z=0)      = %.3f m/s\n', Vmag(is_entrance));
fprintf('Inlet mass flow rate (z=0)     = %.6e kg/s\n', mdot(is_entrance));
fprintf('Outlet mass flow rate          = %.6e kg/s   (%.2f%% of inlet - see ongoing mass-conservation investigation)\n', ...
    mdot(end), 100*mdot(end)/mdot(is_entrance));
fprintf('NOTE: h(z)/Nu(z) also get noisy near the outlet where dT < ~3 K\n');
fprintf('      (small-denominator sensitivity) - check the tail of the table.\n');

results.z       = z;
results.q       = q;
results.Twall   = Twall;
results.Tbulk   = Tbulk;
results.dT      = dT;
results.h       = h;
results.Nu      = Nu;
results.Pbulk   = Pbulk;     % raw CFD gauge static pressure (Fluent's operating-pressure reference)
results.dP_bulk = dP_bulk;   % shifted, inlet-referenced (z=0 -> 0), <= 0 downstream
results.Vmag    = Vmag;      % bulk (mass-weighted average) velocity magnitude
results.mdot    = mdot;      % bulk mass flow rate, sign-corrected to positive
results.is_entrance   = is_entrance;   % applies to h/Nu only - see note above
results.h_bar_simple  = h_bar_simple;
results.Nu_bar_simple = Nu_bar_simple;
results.h_bar_lmtd    = h_bar_lmtd;
results.Nu_bar_lmtd   = Nu_bar_lmtd;

end

%% ---- Local helper: parse a Fluent Surface Integral Report ---- %%
%% (works for any "<surface-name>=<z-value>   <reported-value>" line, ---- %%
%%  regardless of the surface-name prefix used) ---- %%
%%
%% Data rows are located by finding the dashed "----" separator line that
%% marks the end of the header block, rather than by counting header
%% lines - the header is 1 or 2 lines depending on whether the variable
%% name (e.g. "Mass Flow Rate") is short enough to share a line with its
%% units, so a fixed line-count offset isn't reliable across report types.
%% Parsing stops at the SECOND separator, which marks the start of the
%% trailing "Net" summary row.
function [z, val] = read_fluent_report(path)
    lines = splitlines(fileread(path));

    is_separator = ~cellfun(@isempty, regexp(lines, '^[-\s]*-[-\s]*$', 'once'));
    sep_idx = find(is_separator);
    if isempty(sep_idx)
        error('read_fluent_report:noHeaderSeparator', ...
            'Could not find the "----" header separator line in %s.', path);
    end
    header_end = sep_idx(1);
    if numel(sep_idx) >= 2
        data_end = sep_idx(2) - 1;   % stop before the second separator (precedes "Net")
    else
        data_end = numel(lines);     % no second separator found - read to end of file
    end

    z = []; val = [];
    for i = header_end+1 : data_end
        tok = regexp(strtrim(lines{i}), '^\S+=(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)$', 'tokens');
        if ~isempty(tok)
            z(end+1)   = str2double(tok{1}{1}); %#ok<AGROW>
            val(end+1) = str2double(tok{1}{2}); %#ok<AGROW>
        end
    end
    if isempty(z)
        error('read_fluent_report:noDataRows', ...
            'Found the header separator in %s but no data rows after it.', path);
    end

    [z, idx] = sort(z, 'descend');   % z=0 (inlet) -> most negative (outlet)
    val = val(idx);
    z = z(:); val = val(:);
end

function lines = splitlines(text)
    lines = regexp(text, '\r\n|\n', 'split');
end