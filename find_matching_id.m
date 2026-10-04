function id = find_matching_id(inputs)
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

results_dir = 'results';
index_file  = fullfile(results_dir, 'results_index.mat');

id = [];
if ~exist(index_file, 'file')
    return
end

loaded = load(index_file, 'results_index');
results_index = loaded.results_index;

for row = 1:height(results_index)
    if isequal(results_index.inputs{row}, inputs)
        id = results_index.id(row);
        return
    end
end

end