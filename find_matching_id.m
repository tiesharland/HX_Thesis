function id = find_matching_id(inputs, ignore_tol_T)
%FIND_MATCHING_ID  Look up a "physical input configuration" struct in
%results/results_index.mat by exact (isequal) match.
%
%   id = find_matching_id(inputs)
%
% inputs must be built the same way run_disc_model_fwdpass/
% run_lumped_model build their own matching key -- use
% physical_inputs_disc or physical_inputs_lumped to construct it, so a
% caller doing a lookup BEFORE calling the model never builds the key
% differently than the model function itself does when it saves.
%
% Returns the matching id (a double), or [] if results_index.mat doesn't
% exist yet or no row matches. This function only reads the index -- it
% never creates it and never assigns a new id; that happens inside
% run_disc_model_fwdpass/run_lumped_model when they actually save.

% ignore_tol_T (optional, default false): also accept a row that differs only
% in tol_T. Use it for the lumped model, which has no tolerance of its own, so
% it shares the id of the discretised runs of the same physical inputs.
if nargin < 2, ignore_tol_T = false; end

results_dir = 'results';
index_file  = fullfile(results_dir, 'results_index.mat');

id = [];
if ~exist(index_file, 'file')
    return
end

loaded = load(index_file, 'results_index');
results_index = loaded.results_index;

if ~ignore_tol_T
    for row = 1:height(results_index)
        if isequal(results_index.inputs{row}, inputs)
            id = results_index.id(row);
            return
        end
    end
    return
end

% ignore_tol_T: compare with tol_T stripped on both sides, lowest id first,
% so a leftover lumped-only row (e.g. an old id) never wins over the
% discretised row of the same physical inputs.
if ignore_tol_T
    key = strip_tol(inputs);
    for row = 1:height(results_index)
        if isequal(strip_tol(results_index.inputs{row}), key)
            id = results_index.id(row);
            return
        end
    end
end

end

function s = strip_tol(s)
if isfield(s, 'tol_T'), s = rmfield(s, 'tol_T'); end
end