classdef ResultsCollector < handle
    %RESULTSCOLLECTOR  Tiny handle-class accumulator for model evaluations
    %made inside one parfor iteration.
    %
    % Why this exists: a closure like `model_fun = @(N) run_disc_model_fwdpass(...)`
    % can't accumulate results into an ordinary variable in its caller's
    % workspace -- anonymous functions capture by VALUE, not by reference, so
    % each call would just discard its result once find_converged_Nsegments
    % moves on to the next N. A HANDLE object, though, is captured by
    % reference even by an anonymous function, so wrapping each model_fun call
    % as `@(N) capture(run_disc_model_fwdpass(...), collector, ...)` (see the
    % local `capture` function in the sweep scripts that use this) lets every
    % evaluation -- not just the ones find_converged_Nsegments ends up
    % reporting on -- get stashed here as it happens.
    %
    %   collector = ResultsCollector();
    %   collector.add(results, e, r, d2, N_segments, N_cool_seg);
    %   ... (repeat for every model evaluation made in this parfor iteration)
    %   all_entries{i} = collector.entries;   % ordinary sliced parfor output
    %
    % Deliberately does NOT save anything to disk itself -- entries just sit
    % in memory until the parfor loop finishes, at which point the caller
    % should loop over all_entries SERIALLY (on the client, after the parfor)
    % and call save_disc_results on each one. That's what makes this safe to
    % use from inside parfor at all: every worker only ever writes to its own
    % local collector (no shared state between workers), and the only part
    % that touches the shared results/results_index.mat on disk -- the
    % save_disc_results calls -- happens afterward, one at a time, on the
    % client.

    properties
        % Cell array of structs, one per .add() call: results, e, r, d2,
        % N_segments, N_cool_seg -- exactly what save_disc_results needs
        % (together with the physical inputs that stay fixed across an
        % entire sweep, which the caller already has in its own workspace
        % and doesn't need routed through here).
        entries = {}
    end

    methods
        function add(obj, results, e, r, d2, N_segments, N_cool_seg)
            obj.entries{end+1} = struct('results', results, 'e', e, 'r', r, ...
                'd2', d2, 'N_segments', N_segments, 'N_cool_seg', N_cool_seg);
        end
    end

end