%% run_compare_cfd_meshes.m
% Reads the Fluent exports of several mesh levels of the same CFD case with
% process_wall_CFD_data and compares them with the 2D discretised model
% (correlation-based and/or DNS-based, whichever is saved) from the results
% folder, using compare_cfd_meshes_vs_2d.
%
% Folder layout expected for every mesh:
%   <cfd_root>\<mesh folder>\<run folder>\<export file>
% e.g. C:\Users\tiesh\TUDelft\Thesis\CFD\4by4_150ch_wall\run_2024_1\fine\it2000_run\Vmag_air
%
% Meshes whose folder or exports are missing are skipped with a warning, so
% you can list all four now and only the ones that exist get plotted.

clear; clc;

%% ---- Settings ---- %%
cfd_root = 'C:\Users\tiesh\TUDelft\Thesis\CFD\4by4_150ch_wall\run_2024_1';

% One row per mesh, COARSEST first: {legend label, mesh folder, run folder}.
% Change a run folder here once you have named it.
cfd_cases = {
% 'coarse',   'coarse',    'run2';
'baseline', 'baseline',  'solution2/reports';
'fine',     'fine',      'solution2/reports';
'finest',   'fine_v2',   'solution2/reports';
};

% Export files inside each run folder. The first six are required and in the
% order process_wall_CFD_data expects (q, Twall, Tbulk, Pbulk, Vmag, mass flow).
% The rest are optional: density and viscosity (-> Re), wall shear stress
% (-> Cf = |tau|/(0.5 rho V^2)) and y+. Coolant files = the same names with
% 'air' replaced by 'coolant'.
export_files = {'q_air.srp', 'Twall_air.srp', 'Tbulk_air.srp', 'Pbulk_air.srp', 'Vbulk_air.srp', 'm_dot_air.srp'};
opt_names    = {'rho', 'mu', 'tau', 'yplus'};
opt_files    = {'rho_air.srp', 'mu_air.srp', 'tau_wall_air.srp', 'yplus_air.srp'};
include_coolant = true;                  % also compare the coolant side (needs the *_coolant.srp files, plus Tbulk_coolant_in/_out.srp)

% Which saved 2D result(s) to compare with
id         = 4;
N_segments = 50;
N_cool_seg = [0 5 20];             % >= 1
use_DNS    = [false true];     % variants to look for; missing ones are skipped with a warning
include_1d = true;           % also the 1D discretised model (N_cool_seg = 0, same id and N_segments)
include_lumped = true;       % also the lumped model (results/id-<id>/lumped_<id>.mat)

cfd_P_operating = 55000;     % CFD operating pressure [Pa] (Fluent Pbulk is gauge); empty = shift to model inlet pressure

l_index    = 1;              % coolant-side segment compared

%% ---- Hydraulic diameter: read from the saved model results ---- %%
% (equal for air and coolant in the CFD model, so one value serves both)
d_h = [];
res_dir = fullfile('results', sprintf('id-%d', id));
mats = dir(fullfile(res_dir, '*.mat'));
for m = mats'
    ld = load(fullfile(res_dir, m.name), 'results');
    if isfield(ld, 'results') && isfield(ld.results, 'd_h_air')
        d_h = ld.results.d_h_air;
        fprintf('Hydraulic diameter %.15g m, taken from %s\n', d_h, fullfile(res_dir, m.name));
        break
    end
end
if isempty(d_h)
    error('run_compare_cfd_meshes: no saved result with d_h_air found in %s.', res_dir);
end

%% ---- Process the CFD exports (existing meshes only) ---- %%
[cfd, cfd_cool, labels, cool_ends] = load_cfd_exports(cfd_root, cfd_cases, export_files, opt_names, opt_files, d_h, include_coolant);

%% ---- Compare ---- %%
comparison = compare_cfd_meshes_vs_2d(cfd, 'cfd_labels', labels, ...
    'id', id, 'N_segments', N_segments, 'N_cool_seg', N_cool_seg, ...
    'use_DNS', use_DNS, 'include_1d', include_1d, 'include_lumped', include_lumped, 'cfd_coolant', cfd_cool, 'cfd_coolant_ends', cool_ends, 'cfd_P_operating', cfd_P_operating, 'l_index', l_index);