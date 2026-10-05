function out = thermoturb_query_log(mode, entry)
%THERMOTURB_QUERY_LOG  Persistent log of every thermoturb_cached query,
%so you can check afterward how much of a model run's DNS lookups fell
%inside the table's coverage vs. were extrapolated/clamped.
%
%   thermoturb_query_log('log', entry)   -- append one query. Called by
%                                            thermoturb_cached itself;
%                                            you shouldn't need to call
%                                            this directly.
%   thermoturb_query_log('new_set')      -- start a new "set": every
%                                            subsequent 'log' call is
%                                            tagged with the next set id,
%                                            until the next 'new_set'.
%                                            Called once per Q_pred_L_disc
%                                            evaluation by HX_design1_disc
%                                            (one HX length/Q evaluation --
%                                            the natural "one design run"
%                                            grouping, whether that's a
%                                            single 2D direct-marching
%                                            pass over N_segments*N_cool_seg
%                                            cells, or a 1D pass whose
%                                            per-segment inner convergence
%                                            loop makes a variable number
%                                            of queries). You shouldn't
%                                            need to call this directly.
%   thermoturb_query_log('new_iter')     -- start a new, coarser "iter":
%                                            every subsequent 'log' call
%                                            is tagged with the next iter
%                                            id, until the next
%                                            'new_iter'. Called once per
%                                            forward_pass evaluation by
%                                            run_disc_model_fwdpass (one
%                                            mass-flow/diff_P outer
%                                            iteration of the pressure
%                                            balance loop -- each one
%                                            contains however many 'set's
%                                            HX_design1_disc's length
%                                            search made). You shouldn't
%                                            need to call this directly.
%   thermoturb_query_log('new_seg')      -- start a new "seg": every
%                                            subsequent 'log' call is
%                                            tagged with the next seg id,
%                                            until the next 'new_seg'.
%                                            Called once per air segment
%                                            (the shared `for k =
%                                            1:N_segments` loop in
%                                            HX_design1_disc, for BOTH the
%                                            1D and 2D paths) -- one seg
%                                            covers every query made while
%                                            computing that segment,
%                                            including all of the 1D
%                                            path's inner convergence-loop
%                                            passes or all of the 2D
%                                            path's N_seg_cool columns.
%                                            You shouldn't need to call
%                                            this directly.
%   thermoturb_query_log('new_inner')    -- start a new, finer "inner":
%                                            every subsequent 'log' call
%                                            is tagged with the next inner
%                                            id, until the next
%                                            'new_inner'. Called once per
%                                            pass of the 1D path's
%                                            per-segment inner convergence
%                                            loop (calculate_segment_1d's
%                                            `for inner = 1:max_cool`) --
%                                            the 2D path never calls this,
%                                            so its rows all keep whatever
%                                            inner id was current when the
%                                            run started (making every 2D
%                                            row trivially its own
%                                            segment's "last" inner pass,
%                                            since there's only ever one).
%                                            You shouldn't need to call
%                                            this directly.
%   thermoturb_query_log('new_col')      -- start a new, finer "col": every
%                                            subsequent 'log' call is
%                                            tagged with the next col id,
%                                            until the next 'new_col'.
%                                            Called once per coolant column
%                                            of the 2D path's per-segment
%                                            direct march (calculate_segment_2d's
%                                            `for l = 1:N_seg_cool`) -- the
%                                            1D path never calls this, so
%                                            its rows all keep whatever col
%                                            id was current when the run
%                                            started (making every 1D row
%                                            trivially its own segment's
%                                            only column). You shouldn't
%                                            need to call this directly.
%   log = thermoturb_query_log('get')    -- return everything logged so
%                                            far, as a table with columns
%                                            Re, Pr, Tb, Tw, regime,
%                                            Re_extrap, Pr_extrap,
%                                            Tb_extrap, Tw_extrap,
%                                            any_extrap, set, iter, seg,
%                                            inner, col. set is a 1-based
%                                            integer identifying which
%                                            Q_pred_L_disc call each query
%                                            came from (one "design run");
%                                            iter is a coarser 1-based
%                                            integer identifying which
%                                            mass-flow/diff_P outer
%                                            iteration it came from --
%                                            every set belongs to exactly
%                                            one iter. seg identifies
%                                            which air segment a query
%                                            came from (every seg belongs
%                                            to exactly one set); inner
%                                            identifies which pass of that
%                                            segment's 1D inner
%                                            convergence loop it came from
%                                            (constant, and therefore a
%                                            no-op to filter on, for the
%                                            2D path's queries); col
%                                            identifies which coolant
%                                            column of that segment's 2D
%                                            direct march it came from
%                                            (constant, and therefore a
%                                            no-op to filter on, for the
%                                            1D path's queries).
%   thermoturb_query_log('reset')        -- clear the log and reset the
%                                            set/iter/seg/inner/col
%                                            counters. Call this before a
%                                            run if you only want that
%                                            run's queries (otherwise it
%                                            accumulates across every
%                                            thermoturb_cached call made
%                                            this session).
%
% entry (for 'log'): a 1x1 struct with scalar fields Re, Pr, Tb, Tw,
% regime (0 = laminar, 1 = turbulent, 2 = transition blend of both), and
% the four logical *_extrap fields -- true means that dimension fell
% outside the DNS table's coverage for whichever regime(s) were actually
% evaluated for this query. The set/iter/seg/inner/col ids are NOT part of
% entry -- they're stamped automatically from this function's own
% persistent counters, which
% 'new_set'/'new_iter'/'new_seg'/'new_inner'/'new_col' advance.
%
% Storage is preallocated and doubled on growth (not a table that grows
% row-by-row), since thermoturb_cached may call this on every segment of
% every model run.

persistent Re_log Pr_log Tb_log Tw_log regime_log ...
           Re_extrap_log Pr_extrap_log Tb_extrap_log Tw_extrap_log ...
           set_log set_id iter_log iter_id seg_log seg_id inner_log inner_id ...
           col_log col_id n cap

if isempty(n)
    n       = 0;
    cap     = 1024;
    % Start at 0, not 1: 'new_set'/'new_iter'/'new_seg'/'new_inner'/
    % 'new_col' all advance BEFORE the first query of the run is logged,
    % so the first real id of each must come out as 1, not 2 -- starting
    % at 1 here would skip id 1 entirely.
    set_id   = 0;
    iter_id  = 0;
    seg_id   = 0;
    inner_id = 0;
    col_id   = 0;
    Re_log        = zeros(cap, 1);
    Pr_log        = zeros(cap, 1);
    Tb_log        = zeros(cap, 1);
    Tw_log        = zeros(cap, 1);
    regime_log    = zeros(cap, 1);
    Re_extrap_log = false(cap, 1);
    Pr_extrap_log = false(cap, 1);
    Tb_extrap_log = false(cap, 1);
    Tw_extrap_log = false(cap, 1);
    set_log       = zeros(cap, 1);
    iter_log      = zeros(cap, 1);
    seg_log       = zeros(cap, 1);
    inner_log     = zeros(cap, 1);
    col_log       = zeros(cap, 1);
end

switch mode
    case 'new_set'
        set_id = set_id + 1;
        out    = [];

    case 'new_iter'
        iter_id = iter_id + 1;
        out     = [];

    case 'new_seg'
        seg_id = seg_id + 1;
        out    = [];

    case 'new_inner'
        inner_id = inner_id + 1;
        out      = [];

    case 'new_col'
        col_id = col_id + 1;
        out    = [];

    case 'log'
        n = n + 1;
        if n > cap
            cap = cap * 2;
            Re_log(end+1:cap, 1)        = 0;
            Pr_log(end+1:cap, 1)        = 0;
            Tb_log(end+1:cap, 1)        = 0;
            Tw_log(end+1:cap, 1)        = 0;
            regime_log(end+1:cap, 1)    = 0;
            Re_extrap_log(end+1:cap, 1) = false;
            Pr_extrap_log(end+1:cap, 1) = false;
            Tb_extrap_log(end+1:cap, 1) = false;
            Tw_extrap_log(end+1:cap, 1) = false;
            set_log(end+1:cap, 1)       = 0;
            iter_log(end+1:cap, 1)      = 0;
            seg_log(end+1:cap, 1)       = 0;
            inner_log(end+1:cap, 1)     = 0;
            col_log(end+1:cap, 1)       = 0;
        end
        Re_log(n)        = entry.Re;
        Pr_log(n)        = entry.Pr;
        Tb_log(n)        = entry.Tb;
        Tw_log(n)        = entry.Tw;
        regime_log(n)    = entry.regime;
        Re_extrap_log(n) = entry.Re_extrap;
        Pr_extrap_log(n) = entry.Pr_extrap;
        Tb_extrap_log(n) = entry.Tb_extrap;
        Tw_extrap_log(n) = entry.Tw_extrap;
        set_log(n)       = set_id;
        iter_log(n)      = iter_id;
        seg_log(n)       = seg_id;
        inner_log(n)     = inner_id;
        col_log(n)       = col_id;
        out = [];

    case 'get'
        Re        = Re_log(1:n);
        Pr        = Pr_log(1:n);
        Tb        = Tb_log(1:n);
        Tw        = Tw_log(1:n);
        regime    = regime_log(1:n);
        Re_extrap = Re_extrap_log(1:n);
        Pr_extrap = Pr_extrap_log(1:n);
        Tb_extrap = Tb_extrap_log(1:n);
        Tw_extrap = Tw_extrap_log(1:n);
        set       = set_log(1:n);
        iter      = iter_log(1:n);
        seg       = seg_log(1:n);
        inner     = inner_log(1:n);
        col       = col_log(1:n);
        any_extrap = Re_extrap | Pr_extrap | Tb_extrap | Tw_extrap;
        out = table(Re, Pr, Tb, Tw, regime, Re_extrap, Pr_extrap, Tb_extrap, Tw_extrap, any_extrap, set, iter, seg, inner, col);

    case 'reset'
        n        = 0;
        set_id   = 0;   % see the isempty(n) block above for why 0, not 1
        iter_id  = 0;
        seg_id   = 0;
        inner_id = 0;
        col_id   = 0;
        out      = [];

    otherwise
        error('thermoturb_query_log:UnknownMode', 'Unknown mode "%s".', mode);
end

end