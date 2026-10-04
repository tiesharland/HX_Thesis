function fname = disc_results_filename(N_segments, N_cool_seg, use_DNS, id, solve_T)
%DISC_RESULTS_FILENAME  Results filename convention for the discretised
%HX model (run_disc_model_fwdpass).
%
%   N_cool_seg == 0  -> a<N_segments>_<id>               (old iterative 1D)
%   N_cool_seg > 0   -> a<N_segments>-c<N_cool_seg>_<id>  (2D discretised)
%
% "a" prefixes the air-side segment count (N_segments), "c" prefixes the
% coolant-side segment count (N_cool_seg). "_dns" is inserted right
% before "_<id>" when use_DNS is true, and dropped entirely (no
% placeholder) when false.
%
% solve_T (optional, default false): appends "-T" right after "_<id>"
% when true, flagging that this run solved for the exit temperature
% rather than for the specified Q_tot. Any version suffix from
% next_versioned_filename (e.g. "-1") is applied afterward, so it lands
% after "-T", e.g. "a50_3-T-1.mat".
%
% Shared between run_disc_model_fwdpass (when it saves) and any calling
% script checking, before calling it, whether a given model
% variant's results file already exists for a given id.

if nargin < 5
    solve_T = false;
end

if N_cool_seg == 0
    base = sprintf('a%d', N_segments);
else
    base = sprintf('a%d-c%d', N_segments, N_cool_seg);
end
if use_DNS
    base = [base '_dns'];
end
fname = sprintf('%s_%d.mat', base, id);

if solve_T
    [~, name, ext] = fileparts(fname);
    fname = sprintf('%s-T%s', name, ext);
end

end