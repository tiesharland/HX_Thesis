function fname = lumped_results_filename(id, solve_T)
%LUMPED_RESULTS_FILENAME  Results filename convention for the lumped HX
%model (run_lumped_model).
%
%   lumped_<id>.mat          (solve_T false/omitted)
%   lumped_<id>-T.mat        (solve_T true)
%
% solve_T (optional, default false): appends "-T" right after "_<id>"
% when true, flagging that this run solved for the exit temperature
% rather than for the specified Q_tot. Any version suffix from
% next_versioned_filename (e.g. "-1") is applied afterward, so it lands
% after "-T", e.g. "lumped_3-T-1.mat".
%
% Shared between run_lumped_model (when it saves) and any calling script
% checking, before calling it, whether a results file already exists for
% a given id.

if nargin < 2
    solve_T = false;
end

fname = sprintf('lumped_%d.mat', id);

if solve_T
    [~, name, ext] = fileparts(fname);
    fname = sprintf('%s-T%s', name, ext);
end

end