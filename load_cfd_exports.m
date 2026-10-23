function [cfd, cfd_cool, labels, cool_ends] = load_cfd_exports(cfd_root, cfd_cases, export_files, opt_names, opt_files, d_h, include_coolant, cum_file)
% load_cfd_exports  Read the Fluent exports of several CFD meshes with
% process_wall_CFD_data (existing meshes only).
%
%   cfd_root     : folder holding one subfolder per mesh
%   cfd_cases    : n x 3 cell, COARSEST first: {legend label, mesh folder, run folder}
%   export_files : the six required air-side files (q, Twall, Tbulk, Pbulk, Vmag, mass flow)
%   opt_names / opt_files : optional exports as name-value pairs for
%                  process_wall_CFD_data, e.g. {'rho','mu','tau','yplus'} and
%                  {'rho_air.srp','mu_air.srp','tau_wall_air.srp','yplus_air.srp'}.
%                  Missing optional files are left out with a warning.
%   d_h          : hydraulic diameter [m] (same for air and coolant here)
%   include_coolant : also read the coolant files (the same names with 'air'
%                  replaced by 'coolant'); switched off for all meshes if any
%                  mesh lacks them.
%   cum_file     : (optional, default 'cum_heat_transfer.txt') file written by the Fluent UDF
%                  cum_heat_transfer_udf.c, in each run folder next to the other exports
%                  (columns z, slice heat, cumulative heat; '%' lines are comments). If present,
%                  it is stored as cfd{k}.Q_cum_wall (.z, .Q) and used for the CFD cumulative heat
%                  transfer; otherwise the energy balance on the air is used.
%
% Returns cfd (air), cfd_cool (coolant, {} if off) and labels of the meshes read.
% A mesh missing a required air export is skipped with a warning; if none
% can be read the function stops with an error.
%
% cool_ends (4th output, {} if the coolant is off): per mesh a struct with the coolant
% bulk temperature on the coolant INLET and OUTLET surfaces, read from the files
% '<Tbulk name>_in' / '<Tbulk name>_out' (e.g. Tbulk_coolant_in.srp, Tbulk_coolant_out.srp),
% each at its own stations:
%     z_in, Tin, z_out, Tout       (always, else the struct is [])
%     z_dP, dP                     (only if Pbulk_coolant_in / _out exist: P_in - P_out)
% The complete-flow coolant files (cfd_cool) are the normal export names with
% 'air' replaced by 'coolant'. Wall quantities (heat flux, shear, y+) do not exist on the
% inlet / outlet surfaces, so only the bulk values are read there.

if nargin < 8 || isempty(cum_file), cum_file = 'cum_heat_transfer.txt'; end
cool_files     = strrep(export_files, 'air', 'coolant');
cool_opt_files = strrep(opt_files,   'air', 'coolant');
cfd = {};  cfd_cool = {};  labels = {};  cool_ok = true;  cool_ends = {};

for k = 1:size(cfd_cases, 1)
    run_dir = fullfile(cfd_root, cfd_cases{k,2}, cfd_cases{k,3});
    paths   = fullfile(run_dir, export_files);
    have    = cellfun(@isfile, paths);
    if ~all(have)
        warning('Mesh "%s" skipped -- missing required export(s) in %s:\n  %s', ...
            cfd_cases{k,1}, run_dir, strjoin(export_files(~have), ', '));
        continue
    end
    args = optional_args(run_dir, opt_names, opt_files, cfd_cases{k,1});
    fprintf('\n===== CFD case: %s  (%s) =====\n', cfd_cases{k,1}, run_dir);
    cfd{end+1} = process_wall_CFD_data(paths{:}, d_h, args{:}); %#ok<AGROW>
    labels{end+1} = cfd_cases{k,1};                            %#ok<AGROW>
    cum_path = fullfile(run_dir, cum_file);
    if isfile(cum_path)
        [cz, cQ] = read_cum_heat(cum_path);
        cfd{end}.Q_cum_wall = struct('z', cz, 'Q', cQ);
        Q_bal = cumulative_heat_transfer(cfd{end}, -cz(end));
        fprintf('Cumulative heat (UDF wall flux, %s): %.4f W at z = %.3f m; energy balance on the air: %.4f W\n', ...
            cum_file, cQ(end), cz(end), Q_bal);
    else
        fprintf('No %s in this run folder -- cumulative heat from the energy balance on the air.\n', cum_file);
    end
    if include_coolant
        cpaths = fullfile(run_dir, cool_files);
        chave  = cellfun(@isfile, cpaths);
        if ~all(chave)
            warning('Mesh "%s": coolant export(s) missing (%s) -- coolant comparison switched off.', ...
                cfd_cases{k,1}, strjoin(cool_files(~chave), ', '));
            cool_ok = false;
        else
            cargs = optional_args(run_dir, opt_names, cool_opt_files, cfd_cases{k,1});
            fprintf('--- coolant side ---\n');
            cfd_cool{end+1} = process_wall_CFD_data(cpaths{:}, d_h, cargs{:}, 'k_fun', @k_EG); %#ok<AGROW>
            cool_ends{end+1} = read_cool_ends(run_dir, cool_files{3}, cool_files{4}, cfd_cases{k,1}); %#ok<AGROW>
        end
    end
end
if ~cool_ok || numel(cfd_cool) ~= numel(cfd), cfd_cool = {}; cool_ends = {}; end
if ~isempty(cool_ends) && all(cellfun(@isempty, cool_ends)), cool_ends = {}; end
if isempty(cfd)
    error('load_cfd_exports: none of the listed meshes has a complete set of exports under %s.', cfd_root);
end
fprintf('\nCFD meshes read: %s\n', strjoin(labels, ', '));
end

function e = read_cool_ends(run_dir, T_file, P_file, label)
% Coolant inlet / outlet surface values: '<name>_in' and '<name>_out' next to the others.
e = [];
[T_in_f, T_out_f] = end_names(run_dir, T_file);
if ~isfile(T_in_f) || ~isfile(T_out_f)
    warning('Mesh "%s": coolant inlet/outlet temperature exports not found (%s, %s) -- inlet/outlet figures left out.', ...
        label, T_in_f, T_out_f);
    return
end
[e.z_in,  e.Tin]  = read_fluent_report(T_in_f);
[e.z_out, e.Tout] = read_fluent_report(T_out_f);
[P_in_f, P_out_f] = end_names(run_dir, P_file);
if isfile(P_in_f) && isfile(P_out_f)
    [zi, Pi] = read_fluent_report(P_in_f);
    [zo, Po] = read_fluent_report(P_out_f);
    [e.z_dP, ia, ib] = intersect(zi, zo, 'stable');
    e.dP = Pi(ia) - Po(ib);
else
    warning('Mesh "%s": coolant inlet/outlet pressure exports not found -- coolant pressure-drop figure left out.', label);
end
end

function [f_in, f_out] = end_names(run_dir, name)
[~, n, ext] = fileparts(name);
f_in  = fullfile(run_dir, [n '_in'  ext]);
f_out = fullfile(run_dir, [n '_out' ext]);
end

function args = optional_args(run_dir, names, files, label)
args = {};
for i = 1:numel(names)
    f = fullfile(run_dir, files{i});
    if isfile(f)
        args = [args, {names{i}, f}]; %#ok<AGROW>
    else
        warning('Mesh "%s": optional export %s not found -- left out.', label, files{i});
    end
end
end

function [z, val] = read_fluent_report(path)
% Parse a Fluent Surface Integral Report (.srp): '<surface>=<z>  <value>' rows between the
% first and second dashed separator. Same parser as in process_wall_CFD_data.m.
% Returns column vectors, z sorted descending (z = 0 = inlet first).

lines = regexp(fileread(path), '\r\n|\n', 'split');

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

[z, idx] = sort(z, 'descend');   % z = 0 (inlet) -> most negative (outlet)
val = val(idx);
z = z(:); val = val(:);
end

function [z, Q] = read_cum_heat(path)
% Fluent UDF output: lines starting with '%' are comments, data rows are z, slice [W], cumulative [W].
lines = regexp(fileread(path), '\r\n|\n', 'split');
lines = strtrim(lines);
lines = lines(~cellfun(@isempty, lines) & ~startsWith(lines, '%'));
v = sscanf(strjoin(lines, ' '), '%f', [3 Inf])';
z = v(:,1);  Q = v(:,3);
end