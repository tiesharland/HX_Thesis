function results_file = save_disc_results(results, e, r, hx_theta, fan, fpr_init, ...
    M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
    d2_init, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, ...
    solve_T, tol_T, N_segments, N_cool_seg, use_DNS, overwrite)
%SAVE_DISC_RESULTS  Save one run_disc_model_fwdpass results struct to disk
%under the shared results/results_index.mat convention.
%
% Split out of run_disc_model_fwdpass itself so a caller that already has
% a COMPUTED results struct in hand -- e.g. one collected from a parfor
% worker via a ResultsCollector, rather than one it just computed inline
% -- can save it through the exact same id-matching/versioned-filename
% logic, without re-running the model just to get it saved.
% run_disc_model_fwdpass's own 'save_results' branch now just calls this.
%
% Every argument except `results`, `N_segments`, `N_cool_seg`, `use_DNS`
% and `overwrite` is exactly the physical input you'd pass into
% run_disc_model_fwdpass to (re)compute this same results struct --
% physical_inputs() below is built from them to find/assign this run's
% shared "physical input configuration" id, same as run_disc_model_fwdpass
% does internally when it saves.
%
% Results files live under results/id-<id>/ -- one subfolder per physical
% input configuration -- rather than flat in results/ itself, so a
% sweep's many resolution variants for the same design point sit
% together: results/id-7/a50_7.mat, results/id-7/a50-c5_7-1.mat, etc.
% results_index.mat itself stays directly under results/ (it has to: it's
% the single lookup shared across every id, not something that belongs to
% one of them).
%
% overwrite: true  -> overwrite that id+variant's results file in place.
%            false -> keep existing file(s), save this run as a new
%            version alongside them (a50_3.mat, then a50_3-1.mat,
%            a50_3-2.mat, ...). A name only actually gets a "-N" suffix
%            when that exact file is already there (next_versioned_filename
%            returns the plain name untouched otherwise) -- so a single
%            clean run, saving each (id, N_segments, N_cool_seg) variant
%            at most once, never produces "-1"-suffixed files; those only
%            start appearing once you re-run something that's already
%            been saved before.
%
% NOT safe to call concurrently from multiple parfor/parpool workers -- it
% does an unlocked read-modify-write of the shared results_index.mat (same
% as run_disc_model_fwdpass's own save block always has). Call it only
% from the client/main process: serially, e.g. in a plain `for` loop after
% a parfor sweep has finished computing (and collected) everything it
% needs saved -- never from inside a parfor loop body itself.
%
% Returns the path actually saved to.

results_dir = 'results';
index_file  = fullfile(results_dir, 'results_index.mat');

if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% Physical input configuration shared across ALL model variants
% (everything fwdpass takes EXCEPT the model-selector trio
% N_segments/N_cool_seg/use_DNS, plus solve_T) -- built by the shared
% helper so a caller doing its own cache lookup (via find_matching_id)
% before calling this function builds the exact same key.
inputs = physical_inputs(e, r, hx_theta, fan, fpr_init, ...
    M_dot_coolant, n_modules, Q_tot, T_in_fc, T_out_fc, h, p11, t11, ...
    d2_init, AR_diff, AR_noz, V_inf, R, flight_phase, M_dot_FOD, M_dot_comp, ...
    solve_T, tol_T);

if exist(index_file, 'file')
    loaded = load(index_file, 'results_index');
    results_index = loaded.results_index;
else
    results_index = table('Size', [0 2], 'VariableTypes', {'double', 'cell'}, ...
        'VariableNames', {'id', 'inputs'});
end

% Match against existing entries so the same physical-input case always
% reuses the same id (and therefore the same results filename), across
% however many times it gets saved.
id = [];
for row = 1:height(results_index)
    if isequal(results_index.inputs{row}, inputs)
        id = results_index.id(row);
        break
    end
end

if isempty(id)
    id = height(results_index) + 1;
    new_row = table(id, {inputs}, 'VariableNames', {'id', 'inputs'});
    results_index = [results_index; new_row]; %#ok<AGROW>
    save(index_file, 'results_index');
end

id_dir = fullfile(results_dir, sprintf('id-%d', id));
if ~exist(id_dir, 'dir')
    mkdir(id_dir);
end

base_fname = disc_results_filename(N_segments, N_cool_seg, use_DNS, id, solve_T);

if overwrite
    results_file = fullfile(id_dir, base_fname);
else
    results_file = fullfile(id_dir, ...
        next_versioned_filename(id_dir, base_fname));
end

save(results_file, 'results');
fprintf('Saved results: %s\n', results_file);

end