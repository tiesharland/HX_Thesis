function fname = next_versioned_filename(results_dir, base_fname)
%NEXT_VERSIONED_FILENAME  First free "-<n>" suffixed filename for
%base_fname inside results_dir (used when a model function's hard-coded
%overwrite flag is false).
%
%   fname = next_versioned_filename(results_dir, base_fname)
%
% Returns base_fname unchanged if it doesn't exist yet (first save of
% this case). Otherwise tries "-1", "-2", ... suffixes (inserted before
% the extension) until a free one is found, e.g. "a50_3.mat" ->
% "a50_3-1.mat" -> "a50_3-2.mat".
%
% Shared between run_disc_model_fwdpass and run_lumped_model.

if ~exist(fullfile(results_dir, base_fname), 'file')
    fname = base_fname;
    return
end

[~, name, ext] = fileparts(base_fname);
v = 1;
while exist(fullfile(results_dir, sprintf('%s-%d%s', name, v, ext)), 'file')
    v = v + 1;
end
fname = sprintf('%s-%d%s', name, v, ext);

end