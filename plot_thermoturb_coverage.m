function plot_thermoturb_coverage(varargin)
%PLOT_THERMOTURB_COVERAGE  Visualize how a model run's DNS lookups sit
%relative to the thermoturb table's coverage.
%
%   plot_thermoturb_coverage()                 -- 1 panel, 3D scatter + wireframe box (default)
%   plot_thermoturb_coverage('mode','2d')      -- 3 panels of 2D projections
%   plot_thermoturb_coverage('highlight', K)   -- color the LAST K logged
%                                                  queries bright red
%                                                  (e.g. K = N_segments*N_seg_cool
%                                                  for the final design
%                                                  iteration of a sweep),
%                                                  drawn on top of the rest.
%                                                  K = 0 (default): off,
%                                                  no query stands out.
%                                                  Either way the bulk
%                                                  points are colored by
%                                                  which design run they
%                                                  belong to (see below),
%                                                  on MATLAB's default
%                                                  colormap.
%   plot_thermoturb_coverage('log', log_tbl)   -- plot from a previously
%                                                  saved query log (e.g.
%                                                  results.thermoturb_log
%                                                  from a saved results
%                                                  file) instead of the
%                                                  live persistent log.
%                                                  Lets you replot any old
%                                                  run's coverage without
%                                                  re-running the model.
%   plot_thermoturb_coverage('iter', M)        -- start the plot already
%                                                  filtered to mass-flow
%                                                  iteration M (1 = first,
%                                                  increasing from there --
%                                                  this argument is still
%                                                  the ABSOLUTE iteration
%                                                  number; the dropdown it
%                                                  feeds into uses the
%                                                  offset scheme below).
%                                                  Default [] shows all
%                                                  iterations.
%
%   All three filter dropdowns on the figure -- "Mass-flow iteration",
%   "Q_pred,L call #" (appears whenever the log has both 'set' and 'iter'
%   columns), and "Inner iterations" (appears whenever it has both 'seg'
%   and 'inner') -- share the same scheme: 'All', 'Last (converged)'
%   (offset 0: that level's own final/converged value), then '-1','-2',
%   ... counted BACKWARD from there (one before the last, two before,
%   ...). Picking '-2' in a dropdown always means "two before whatever
%   this level settled on", lining up the SAME relative position across
%   every group at that level regardless of how long each one individually
%   took to get there -- e.g. "Inner iterations" = '-2' with "Q_pred,L
%   call #" = 'All' shows the 3rd-to-last pass of every air segment's
%   convergence loop across every length-search call at once, and a
%   segment that happened to need one extra pass is the only one that
%   reaches one option further back than the rest. The three combine
%   (iter AND call AND inner): e.g. iteration 12 + call 'Last (converged)'
%   + inner '-3' isolates exactly the 4th-to-last pass of iteration 12's
%   final length-search trial. For 2D-path logs (no inner convergence loop
%   at all) every "Inner iterations" choice is identical, since every
%   segment only ever made one pass anyway.
%
%   The call-# and inner dropdowns' NUMBER options aren't fixed -- they're
%   rebuilt after every change to only what's actually reachable from the
%   current coarser selection(s) (e.g. once you've picked one iteration,
%   the call-# list only offers the offsets THAT iteration's length search
%   reached; once you've also picked one call, the inner-pass list only
%   offers the offsets THAT one design run's segments reached). A
%   selection that's still valid after a coarser dropdown changes is kept;
%   one that no longer exists falls back to that dropdown's default
%   ('All' for iter/call, 'Last (converged)' for inner). A "N / total
%   points shown" counter next to the dropdowns updates with any
%   selection.
%
%   Two more dropdowns, on a second row, filter by physical position
%   rather than by convergence loop, so they use plain forward numbering
%   ('All', '1', '2', ...) instead of the backward-from-converged scheme
%   above:
%     - "Air segment" -- which air segment (1..N_segments) a query came
%       from. Shown for BOTH 1D- and 2D-path logs, since both call
%       thermoturb_query_log('new_seg') once per air segment identically.
%       Picking e.g. '3' shows every query from the 3rd air segment,
%       across whatever mass-flow iteration/call/coolant-segment is also
%       selected.
%     - "Coolant segment" -- which coolant column (1..N_seg_cool) a query
%       came from. Only meaningful for a 2D-path log (calculate_segment_2d
%       calls thermoturb_query_log('new_col') once per column); a 1D-path
%       log never advances this, so the dropdown doesn't appear at all for
%       one. Picking e.g. '2' shows every query from the 2nd coolant
%       column, across whatever air segments are also selected.
%   "Inner iterations" (1D-path only) and "Coolant segment" (2D-path only)
%   therefore never both appear at once -- whichever path produced the log
%   determines which of the two shows up, automatically, alongside "Air
%   segment" which always shows (when the log has a 'seg' column at all).
%   Air segment and Coolant segment are independent of each other (they're
%   two different axes of the same k-by-l grid, not one nested inside the
%   other), but both are still rescoped by the coarser iter/call/inner
%   selections the same way call/inner are.
%
%   Run this after a model run (anything that calls thermoturb_cached
%   with log_queries = true). By default it reads the accumulated log
%   from thermoturb_query_log('get'); pass 'log' to use a saved one
%   instead. Either way, the table's own coverage bounds are read
%   straight from thermoturb_table.mat.
%
%   Coloring: points are colored by "design run" rather than by
%   individual query. Every thermoturb_cached call made during one
%   Q_pred_L_disc evaluation (one HX length/Q evaluation -- HX_design1_disc
%   marks this via thermoturb_query_log('new_set')) shares one color,
%   whether that evaluation made exactly N_segments*N_cool_seg queries
%   (2D direct-marching, one per segment-column cell) or a variable
%   number (1D, whose per-segment inner convergence loop can take a
%   different number of iterations each run -- the resulting queries are
%   "progression of the average over time", not a physical progression
%   through the coolant channel, but they still all belong to the same
%   Q_pred_L output and so get the same color).
%
%   The color VALUE is each set's chronological rank WITHIN its own
%   mass-flow iteration (1st length-search call of that iteration, 2nd,
%   ...), not a raw running count over the whole run -- so looking at one
%   mass-flow iteration shows exactly how many Q_pred_L_disc calls it
%   took to find L_solution, as a clean small gradient; looking at
%   several/all iterations together shows that same small gradient
%   repeat once per iteration. The color axis is rescaled to the current
%   selection every time (initial draw and every dropdown change), so
%   the colormap always stretches to fit whatever's on screen rather than
%   staying fixed to the whole run's range. Older saved logs (from before
%   the 'iter' column existed) fall back to a plain rank of 'set' over
%   the whole log; logs from before 'set' existed fall back further to
%   per-query chronological coloring, both with a warning.
%
%   "Coverage" here means: a query's (Pr,Tb,Tw) cell is inside the
%   sampled grid. Re is NOT part of that test -- within an in-grid
%   (Pr,Tb,Tw) combo, the whole laminar -> transition-blend ->
%   turbulent-extrapolated -> turbulent-DNS path is considered covered,
%   since thermoturb_cached always returns a value there (it's only
%   extrapolating in Re, same combo). Only a (Pr,Tb,Tw) that falls
%   outside the sampled grid entirely is "not covered".
%
%   Coverage box: not every (Pr,Tb,Tw) combo's Re sweep spans the same
%   range, and Pr varies the range only slightly, so rather than drawing
%   one flat min/max rectangle, the box is built from the corners of the
%   sampled (Tb,Tw) footprint -- the vertices of its convex hull. For a
%   rectangular grid that is the usual 4 corners (min/max Tb x min/max
%   Tw); the table is no longer rectangular in (Tb,Tw) (Tw <= Tb was never
%   queried), so its hull is a different quadrilateral with a slanted edge
%   along Tw = Tb + step. At each hull vertex, across the Pr values
%   sampled there, the Pr-slice with the NARROWEST Re span is used as that
%   vertex's representative range (the conservative choice given Pr's
%   small effect). That gives three Re values per vertex: laminar start,
%   where real turbulent DNS data ends, and turbulent end (incl.
%   extrapolated). Re=2300/3000 (the fixed laminar/turbulent demarcation)
%   are drawn as plain reference lines.
%
%   In 2D, each level is the hull ring projected onto the panel's axes
%   (Re against Tb, or Re against Tw), so a panel shows the boundary of
%   every hull edge at once rather than merging the other temperature away.
%
%   In 3D, the three boundaries form one wireframe box: vertical edges
%   run laminar-start -> DNS-end -> turbulent-end at each hull vertex,
%   with a ring around the hull at each of those three levels. All three
%   levels are wireframe-only -- no filled plane.
%
%   Note: thermoturb_cached flags Tb/Tw extrapolation against the
%   bounding RECTANGLE of the table (min/max Tb, min/max Tw), not this
%   hull, so a query in the empty Tw <= Tb corner of that rectangle is not
%   flagged as Tb_extrap/Tw_extrap even though it sits outside the box
%   drawn here (it gets nearest-neighbour values). Physically it never
%   happens (the wall is hotter than the bulk).
%
%   Mass-flow iteration (iter): one level coarser than the design-run
%   "set" coloring above. run_disc_model_fwdpass marks a new iter via
%   thermoturb_query_log('new_iter') at the top of every forward_pass
%   call -- one mass-flow/diff_P outer iteration of the pressure balance
%   loop, containing however many sets (Q_pred_L_disc calls) the length
%   search inside made at that mass flow. The 'iter' argument and the
%   figure's dropdown filter on this.
%
%   Call thermoturb_query_log('reset') before the run if you only want
%   that run's queries shown (otherwise the log accumulates across every
%   thermoturb_cached call made this session -- including, e.g., a
%   separate Re-sweep test script run earlier).

p = inputParser;
addParameter(p, 'mode', '3d', @(s) any(strcmpi(s, {'2d','3d'})));
addParameter(p, 'highlight', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'log', [], @(t) isempty(t) || istable(t));
addParameter(p, 'iter', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
parse(p, varargin{:});
mode      = lower(p.Results.mode);
highlight = round(p.Results.highlight);
log_input = p.Results.log;
iter_init = p.Results.iter;

% Fixed laminar/turbulent demarcation used by thermoturb_cached -- only
% used here to draw reference lines, not recomputed.
Re_lam_max = 2300;
Re_tur_min = 3000;

%% ---- Load table and compute the coverage boundaries at the hull corners ---- %%

table_file = fullfile(pwd, 'thermoturb_table.mat');
if ~isfile(table_file)
    error('plot_thermoturb_coverage: %s not found.', table_file);
end
loaded = load(table_file, 'tt');
tt     = loaded.tt;

grid_TbTw = unique(tt{:, {'Tb','Tw'}}, 'rows');

% Corners = vertices of the convex hull of the sampled (Tb,Tw) footprint, in
% hull order (consecutive rows are neighbours, so [1:n 1] is the closed
% ring). 'Simplify' drops points that merely lie along an edge. Hull
% vertices are always actual grid points, so every corner has data. For a
% full rectangular grid this is the 4 rectangle corners, as before.
hull_idx = convhull(grid_TbTw(:,1), grid_TbTw(:,2), 'Simplify', true);
corners  = grid_TbTw(hull_idx(1:end-1), :);   % convhull repeats the first point at the end
n_corn   = size(corners, 1);

lam_lo = zeros(n_corn,1);   % laminar start, per corner (box's lower edge)
tur_hi = zeros(n_corn,1);   % turbulent end (incl. extrapolated), per corner (box's upper edge)
tur_lo = zeros(n_corn,1);   % where real turbulent DNS data ends, per corner (middle level)
for i = 1:n_corn
    [lo, ~]    = corner_re_range(tt, corners(i,1), corners(i,2), 0);
    [lo2, hi2] = corner_re_range(tt, corners(i,1), corners(i,2), 1);
    lam_lo(i) = lo;
    tur_lo(i) = lo2;
    tur_hi(i) = hi2;
end

%% ---- Load logged queries ---- %%

if isempty(log_input)
    log_tbl = thermoturb_query_log('get');
else
    log_tbl = log_input;
end
if isempty(log_tbl) || height(log_tbl) == 0
    error(['plot_thermoturb_coverage: the query log is empty. Run the model ' ...
           'first with thermoturb_cached''s log_queries flag set to true, or ' ...
           'pass a saved log via the ''log'' argument (e.g. results.thermoturb_log).']);
end

n_pts = height(log_tbl);

% "Not covered" = the query's (Pr,Tb,Tw) cell fell outside the sampled
% grid. Re_extrap does NOT count -- extrapolating in Re within an
% in-grid combo is covered.
not_covered = log_tbl.Pr_extrap | log_tbl.Tb_extrap | log_tbl.Tw_extrap;

fprintf('plot_thermoturb_coverage: %d logged queries, %d (%.1f%%) outside the (Pr,Tb,Tw) grid.\n', ...
        n_pts, sum(not_covered), 100*mean(not_covered));

% Color by design-run "set" when available: every query made during one
% Q_pred_L_disc call (one HX length/Q evaluation -- the natural "one
% design run", whether it's a 2D pass of N_segments*N_cool_seg cells or a
% 1D pass whose per-segment inner convergence loop makes a variable
% number of queries) shares the same set id and therefore the same
% color, rather than every individual query getting its own shade.
%
% The color VALUE itself is each set's chronological rank WITHIN its own
% mass-flow iteration (1st length-search call, 2nd, ...), not the raw
% globally-increasing set id -- that's what lets one mass-flow iteration
% show "how many evaluations it took to find L_solution" as a clean
% small gradient, the same small gradient repeating across every
% iteration when several/all are shown together. Combined with the
% per-axes caxis('auto') calls below (on first draw and again on every
% dropdown change), the color span always stretches to fit whatever
% subset is actually on screen, so a single filtered-to-one-iteration
% view isn't stuck using only a sliver of a colormap scaled for the
% whole run. Logs saved before 'iter' existed fall back to a plain dense
% rank of 'set' over the whole log; logs from before 'set' existed fall
% back further to per-query chronological order.
if ismember('set', log_tbl.Properties.VariableNames) && ismember('iter', log_tbl.Properties.VariableNames)
    color_idx   = zeros(n_pts, 1);
    call_offset = zeros(n_pts, 1);   % see "Q_pred_L call-number filter" below
    for this_iter = unique(log_tbl.iter)'
        rows = log_tbl.iter == this_iter;
        [~, ~, local_rank] = unique(log_tbl.set(rows));   % ascending -> rank 1,2,3,...
        color_idx(rows)   = local_rank;
        % Counted BACKWARD from that iteration's OWN last call (0 = the
        % call that the length search settled on for this mass-flow
        % iteration), not the forward rank above -- see call_offset's use
        % in the call-number filter for why.
        call_offset(rows) = local_rank - max(local_rank);
    end
    color_label = 'Q_{pred,L} call # within mass-flow iteration';
    % color_idx IS the per-iteration call-rank here, so it doubles as the
    % thing the "Q_pred_L call #" dropdown filters on below.
    color_by_local_rank = true;
elseif ismember('set', log_tbl.Properties.VariableNames)
    [~, ~, color_idx] = unique(log_tbl.set);
    color_label = 'Design run (Q_{pred,L} evaluation)';
    color_by_local_rank = false;
else
    warning(['plot_thermoturb_coverage: this log has no ''set'' column ' ...
             '(saved before set-based coloring was added) -- falling back ' ...
             'to per-query chronological order.']);
    color_idx   = (1:n_pts)';
    color_label = 'Query order (chronological)';
    color_by_local_rank = false;
end

% The trailing 'highlight' queries (e.g. the final design iteration's
% N_segments*N_seg_cool calls) are singled out in bright red and drawn on
% top. The bulk points always stay on MATLAB's default colormap, whether
% or not highlighting is on.
highlight    = min(highlight, n_pts);
is_highlight = false(n_pts, 1);
if highlight > 0
    is_highlight(end-highlight+1:end) = true;
end
bulk_mask     = ~is_highlight;
bulk_cmap     = parula(256);   % MATLAB default, always
col_highlight = [1 0 0];

% Mass-flow iteration filter (coarser than 'set' -- see docstring). Older
% saved logs have no 'iter' column; fall back to showing everything with
% no dropdown, since there's nothing to filter on.
%
% Counted BACKWARD FROM THE LAST mass-flow iteration, same scheme as the
% call-# and inner-iteration filters below: iter_offset = 0 is the run's
% last (converged) mass-flow iteration, -1 the one before that, and so
% on. There's no coarser level above this one, so (unlike call_offset /
% inner_offset) it's a single global offset, not computed per group.
has_iter = ismember('iter', log_tbl.Properties.VariableNames);
if has_iter
    iter_offset = log_tbl.iter - max(log_tbl.iter);
    iter_items  = offset_items_for_scope(iter_offset, true(n_pts, 1));
    if isempty(iter_init)
        iter_sel_idx = 1;   % 'All'
    else
        % iter_init is still an ABSOLUTE iteration number (1 = first) --
        % that's the natural way to ask to start on a specific one -- so
        % convert it to the offset the dropdown actually uses.
        target_offset = iter_init - max(log_tbl.iter);
        match = offset_item_index(iter_items, target_offset);
        if isempty(match)
            warning('plot_thermoturb_coverage: iter=%g not found in the log (max iter = %d) -- showing all.', ...
                     iter_init, max(log_tbl.iter));
            iter_sel_idx = 1;
        else
            iter_sel_idx = match;
        end
    end
    iter_mask_init = offset_mask_from_selection(iter_offset, iter_items{iter_sel_idx});
else
    if ~isempty(iter_init)
        warning(['plot_thermoturb_coverage: this log has no ''iter'' column ' ...
                 '(saved before mass-flow-iteration tracking was added) -- ''iter'' ignored.']);
    end
    iter_mask_init = true(n_pts, 1);
end

% Q_pred_L call-number filter, one level finer than the iter filter above.
% call_offset (computed alongside color_idx above) is each query's call
% rank counted BACKWARD FROM ITS OWN ITERATION'S last call: 0 is the call
% that iteration's length search actually settled on, -1 the one before
% it, and so on -- so filtering on it directly, combined with the iter
% mask, lets you isolate exactly one length-trial's queries within one
% mass-flow iteration (or "every iteration's next-to-last call at once"
% with call_offset == -1 and iter left on 'All', regardless of how many
% total calls each of those iterations individually took). Only offered
% when color_idx/call_offset actually mean that (see color_by_local_rank
% above); for older logs without per-iteration rank there's nothing
% meaningful to filter on, so it's left off entirely.
%
% The dropdown's actual OPTIONS are rebuilt by the callback every time the
% iter dropdown changes, to whatever call offsets really occur within the
% currently-selected iteration(s) -- e.g. iteration 5 might only ever have
% gone through 4 length-search calls (offsets 0..-3), so there's no point
% offering '-4' or beyond just because some other iteration needed more.
% The initial list here uses the same logic against iter_mask_init.
if color_by_local_rank
    call_items     = offset_items_for_scope(call_offset, iter_mask_init);
    call_sel_idx   = 1;   % 'All'
    call_mask_init = true(n_pts, 1);
else
    call_mask_init = true(n_pts, 1);
end

% Inner-iteration filter, finer still: within ONE air segment ('seg'),
% the 1D path's per-segment convergence loop queries thermoturb_cached
% once per pass ('inner') before settling on a converged value -- by
% default only that last, converged pass is worth looking at, with every
% earlier pass just clutter. inner_offset re-numbers 'inner' COUNTING
% BACKWARD FROM CONVERGENCE, same scheme as iter_offset/call_offset above:
% 0 = the converged pass, -1 = one pass before that, -2 = two before, and
% so on (the raw 'inner' id is a single counter running across the whole
% log, so its values mean nothing on their own). This way every segment's
% converged pass lines up on the same number (0) regardless of how many
% passes it took to get there, and a segment that happened to need one
% extra pass is the only one that reaches one number further negative
% than the rest -- exactly the "most segments stop at -3, one outlier
% reaches -4" shape this is meant to show. The 2D path never advances
% 'inner' at all, so every (set,seg) group there has only one row anyway
% -- its offset is always 0 and 'Last'/'All'/any specific number are all
% identical there, which is the correct behaviour (there's no inner loop
% to filter). Only offered when the log actually has 'seg'/'inner'
% columns (both added alongside this feature; older logs fall back to
% showing everything, no dropdown).
%
% Like the call-# dropdown, the actual OPTIONS are rebuilt by the callback
% to whatever occurs within the currently-selected iteration+call (e.g. a
% design run whose segments all converged within 3 passes has no reason
% to offer '-3' or beyond); the initial list here uses the same logic
% against iter_mask_init & call_mask_init.
% Air segment / coolant segment / inner-iteration filters all key off
% 'seg' (and, for the latter two, 'inner'/'col'), so the groupings they
% share are computed once up front.
has_air_seg = ismember('seg', log_tbl.Properties.VariableNames);
if has_air_seg
    group_id_set = findgroups(log_tbl.set);          % one group per set
    group_id_seg = findgroups(log_tbl.set, log_tbl.seg);   % one group per (set,seg)
end

% Inner-iteration filter: within ONE air segment ('seg'), the 1D path's
% per-segment convergence loop queries thermoturb_cached once per pass
% ('inner') before settling on a converged value -- by default only that
% last, converged pass is worth looking at, with every earlier pass just
% clutter. inner_offset re-numbers 'inner' COUNTING BACKWARD FROM
% CONVERGENCE, same scheme as iter_offset/call_offset above: 0 = the
% converged pass, -1 = one pass before that, -2 = two before, and so on
% (the raw 'inner' id is a single counter running across the whole log,
% so its values mean nothing on their own). This way every segment's
% converged pass lines up on the same number (0) regardless of how many
% passes it took to get there, and a segment that happened to need one
% extra pass is the only one that reaches one number further negative
% than the rest -- exactly the "most segments stop at -3, one outlier
% reaches -4" shape this is meant to show.
has_seg_inner = has_air_seg && ismember('inner', log_tbl.Properties.VariableNames);
if has_seg_inner
    max_inner    = splitapply(@max, log_tbl.inner, group_id_seg);
    % 'inner' is assigned by a single global counter that only ever
    % advances (new_inner is called in strict chronological order, with
    % segments themselves processed strictly one at a time), so each
    % group's own ids are a contiguous increasing run -- subtracting off
    % the group's own MAXIMUM turns that straight into a 0-based "distance
    % before convergence" without needing a per-group loop.
    inner_offset = log_tbl.inner - max_inner(group_id_seg);
else
    inner_offset = [];
end
% Only worth a dropdown for a 1D-path run: a 2D-path run never calls
% 'new_inner' at all, so inner_offset is trivially 0 everywhere and every
% choice in the dropdown would behave identically -- just clutter (Air
% segment / Coolant segment below are what matters for 2D instead, and
% 'we don't need to order by inner iteration' there).
show_inner_dropdown = has_seg_inner && any(inner_offset ~= 0);
if show_inner_dropdown
    inner_items     = offset_items_for_scope(inner_offset, iter_mask_init & call_mask_init);
    inner_sel_idx   = 2;   % 'Last (converged)' by default (see offset_items_for_scope)
    inner_mask_init = inner_offset == 0;
else
    inner_mask_init = true(n_pts, 1);
end

% Air segment / coolant segment filters: a PHYSICAL position, not a
% convergence loop, so these use plain forward numbering (1 = first), not
% the backward-from-converged offset scheme above.
%
% Air segment: seg's local rank WITHIN ITS OWN SET, 1..N_segments -- 'seg'
% itself is a single counter running across the whole log (like
% set/iter/inner), so its raw values mean nothing on their own. Offered
% for BOTH the 1D and 2D paths, since both call 'new_seg' once per air
% segment identically: selecting e.g. '3' shows every query from the 3rd
% air segment, across whatever else is currently shown.
if has_air_seg
    min_seg = splitapply(@min, log_tbl.seg, group_id_set);
    % Same contiguous-run argument as inner_offset above, but from the
    % group's MINIMUM (forward rank) rather than its maximum (backward
    % offset) -- segments are processed strictly in order k=1:N_segments
    % within one set, for both paths alike.
    air_seg_rank      = log_tbl.seg - min_seg(group_id_set) + 1;
    air_seg_items     = fwd_items_for_scope(air_seg_rank, iter_mask_init & call_mask_init);
    air_seg_sel_idx   = 1;   % 'All'
    air_seg_mask_init = true(n_pts, 1);
else
    air_seg_mask_init = true(n_pts, 1);
end

% Coolant segment: col's local rank WITHIN ITS OWN (set,seg) GROUP,
% 1..N_seg_cool -- only meaningful for a 2D-path run (calculate_segment_2d
% calls 'new_col' once per coolant column): selecting e.g. '2' shows every
% query from the 2nd coolant column, across whatever air segments are
% currently shown. A 1D-path run never calls 'new_col', so col stays
% constant within every (set,seg) group and the rank collapses to a
% trivial 1 everywhere, in which case there's nothing to offer a dropdown
% for.
has_col = ismember('col', log_tbl.Properties.VariableNames);
if has_col
    min_col = splitapply(@min, log_tbl.col, group_id_seg);
    coolant_seg_rank = log_tbl.col - min_col(group_id_seg) + 1;
else
    coolant_seg_rank = [];
end
show_coolant_dropdown = has_col && any(coolant_seg_rank > 1);
if show_coolant_dropdown
    coolant_seg_items     = fwd_items_for_scope(coolant_seg_rank, iter_mask_init & call_mask_init);
    coolant_seg_sel_idx   = 1;   % 'All'
    coolant_seg_mask_init = true(n_pts, 1);
else
    coolant_seg_mask_init = true(n_pts, 1);
end

%% ---- Colours ---- %%

col_box    = [0.55 0.55 0.55];   % box edges (laminar start / turbulent end)
col_dnsend = [0.55 0.55 0.55];   % turbulent-DNS-ends reference line
col_thresh = [0.55 0.55 0.55];   % Re=2300 / Re=3000 reference lines
col_grid   = [0.55 0.55 0.55];   % Tb-Tw sampled-grid footprint

%% ---- Plot ---- %%

create_highlight = any(is_highlight);
show_bulk_init   = bulk_mask & iter_mask_init & call_mask_init & inner_mask_init & air_seg_mask_init & coolant_seg_mask_init;
show_hl_init     = is_highlight & iter_mask_init & call_mask_init & inner_mask_init & air_seg_mask_init & coolant_seg_mask_init;

if strcmp(mode, '2d')

    fig = figure('Name', 'thermoturb DNS coverage', 'NumberTitle', 'off');

    panels = {
        'Re', 'Tb', 'Re_{bulk}', 'T_b [K]';
        'Re', 'Tw', 'Re_{bulk}', 'T_w [K]';
        'Tb', 'Tw', 'T_b [K]',   'T_w [K]'
    };

    panel_info(3) = struct('ax', [], 'bulk_h', [], 'highlight_h', [], 'x_all', [], 'y_all', [], 'z_all', []);

    for k = 1:3
        ax = subplot(1, 3, k);
        hold(ax, 'on')

        xf = panels{k,1}; yf = panels{k,2};
        x_all = log_tbl.(xf); y_all = log_tbl.(yf);

        % Bulk points first (so reference lines drawn after stay visible
        % on top of dense clusters), colormap set per-axes right after.
        h_bulk = scatter(ax, x_all(show_bulk_init), y_all(show_bulk_init), 12, color_idx(show_bulk_init), 'o', 'filled', ...
                'MarkerFaceAlpha', 0.35, 'HandleVisibility', 'off');
        colormap(ax, bulk_cmap);
        caxis(ax, 'auto');   % rescale color span to whatever's actually shown

        if strcmp(xf, 'Re') && strcmp(yf, 'Tb')
            draw_ring_lines(ax, lam_lo, corners(:,1), col_box,    '-',  1.5, 'Coverage start (laminar)');
            draw_ring_lines(ax, tur_hi, corners(:,1), col_box,    '-',  1.5, 'Coverage end (turbulent, incl. extrap.)');
            draw_ring_lines(ax, tur_lo, corners(:,1), col_dnsend, '--', 1.0, 'Turbulent DNS data ends');
            xline(ax, Re_lam_max, ':', 'Color', col_thresh, 'HandleVisibility', 'off');
            xline(ax, Re_tur_min, ':', 'Color', col_thresh, 'HandleVisibility', 'off');
        elseif strcmp(xf, 'Re') && strcmp(yf, 'Tw')
            draw_ring_lines(ax, lam_lo, corners(:,2), col_box,    '-',  1.5, 'Coverage start (laminar)');
            draw_ring_lines(ax, tur_hi, corners(:,2), col_box,    '-',  1.5, 'Coverage end (turbulent, incl. extrap.)');
            draw_ring_lines(ax, tur_lo, corners(:,2), col_dnsend, '--', 1.0, 'Turbulent DNS data ends');
            xline(ax, Re_lam_max, ':', 'Color', col_thresh, 'HandleVisibility', 'off');
            xline(ax, Re_tur_min, ':', 'Color', col_thresh, 'HandleVisibility', 'off');
        else
            % Sampled (Tb,Tw) footprint = the hull ring itself
            ring = [1:n_corn 1];
            plot(ax, corners(ring,1), corners(ring,2), '--', 'Color', col_grid, 'LineWidth', 1.2, ...
                 'HandleVisibility', 'off');
            plot(ax, nan, nan, '--', 'Color', col_grid, 'LineWidth', 1.2, ...
                 'DisplayName', 'Sampled (Tb,Tw) grid');
        end

        % Highlighted points drawn last, on top of everything.
        h_hl = [];
        if create_highlight
            h_hl = scatter(ax, x_all(show_hl_init), y_all(show_hl_init), 30, col_highlight, 'o', 'filled', ...
                    'MarkerEdgeColor', 'k', 'DisplayName', 'Final design iteration', ...
                    'HandleVisibility', tern(k==1));
        end

        xlabel(ax, panels{k,3}); ylabel(ax, panels{k,4});
        title(ax, sprintf('%s vs %s', panels{k,4}, panels{k,3}))
        grid(ax, 'on')
        cb = colorbar(ax);
        cb.Label.String = color_label;
        if k == 1, legend(ax, 'Location', 'best', 'FontSize', 7); end

        panel_info(k) = struct('ax', ax, 'bulk_h', h_bulk, 'highlight_h', h_hl, ...
                                'x_all', x_all, 'y_all', y_all, 'z_all', []);
    end

    sgtitle(fig, 'DNS query coverage vs. thermoturb table bounds')

else % 3d

    fig = figure('Name', 'thermoturb DNS coverage (3D)', 'NumberTitle', 'off');
    ax = axes(fig);
    hold(ax, 'on')

    h_bulk = scatter3(ax, log_tbl.Re(show_bulk_init), log_tbl.Tb(show_bulk_init), log_tbl.Tw(show_bulk_init), ...
             12, color_idx(show_bulk_init), 'o', 'filled', 'MarkerFaceAlpha', 0.35, 'HandleVisibility', 'off');
    colormap(ax, bulk_cmap);
    caxis(ax, 'auto');   % rescale color span to whatever's actually shown

    draw_wireframe_box(ax, lam_lo, tur_lo, tur_hi, corners, col_box, 'Coverage box (laminar start / DNS end / turbulent end)');
    draw_re_frame(ax, Re_lam_max, corners, col_dnsend, sprintf('Laminar regime ends (Re = %d)', Re_lam_max));
    draw_re_frame(ax, Re_tur_min, corners, col_dnsend, sprintf('Turbulent regime starts (Re = %d)', Re_tur_min));

    h_hl = [];
    if create_highlight
        h_hl = scatter3(ax, log_tbl.Re(show_hl_init), log_tbl.Tb(show_hl_init), log_tbl.Tw(show_hl_init), ...
                 36, col_highlight, 'o', 'filled', 'MarkerEdgeColor', 'k', ...
                 'DisplayName', 'Final design iteration');
    end

    xlabel(ax, 'Re_{bulk}'); ylabel(ax, 'T_b [K]'); zlabel(ax, 'T_w [K]');
    title(ax, 'DNS query coverage vs. thermoturb table bounds')
    cb = colorbar(ax);
    cb.Label.String = color_label;
    legend(ax, 'Location', 'best')
    grid(ax, 'on')
    view(ax, 3)

    panel_info = struct('ax', ax, 'bulk_h', h_bulk, 'highlight_h', h_hl, ...
                         'x_all', log_tbl.Re, 'y_all', log_tbl.Tb, 'z_all', log_tbl.Tw);

end

%% ---- Mass-flow iteration dropdown ---- %%

if has_iter
    info = struct('log_tbl', log_tbl, 'bulk_mask', bulk_mask, 'is_highlight', is_highlight, ...
                   'color_idx', color_idx, 'mode', mode, 'panels', panel_info, ...
                   'iter_offset', iter_offset, 'color_by_local_rank', color_by_local_rank, ...
                   'iter_popup', [], 'call_popup', [], 'inner_popup', [], ...
                   'air_seg_popup', [], 'coolant_seg_popup', [], 'count_text', [], ...
                   'call_offset', [], 'inner_offset', [], 'air_seg_rank', [], 'coolant_seg_rank', []);

    uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
              'Position', [0.01 0.955 0.16 0.03], 'String', 'Mass-flow iteration:', ...
              'HorizontalAlignment', 'left', 'BackgroundColor', fig.Color);
    iter_popup = uicontrol(fig, 'Style', 'popupmenu', 'Units', 'normalized', ...
              'Position', [0.17 0.955 0.1 0.035], 'String', iter_items, ...
              'Value', iter_sel_idx, 'Callback', @thermoturb_coverage_filter_callback);
    info.iter_popup = iter_popup;

    if color_by_local_rank
        uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
                  'Position', [0.28 0.955 0.13 0.03], 'String', 'Q_{pred,L} call #:', ...
                  'HorizontalAlignment', 'left', 'BackgroundColor', fig.Color);
        call_popup = uicontrol(fig, 'Style', 'popupmenu', 'Units', 'normalized', ...
                  'Position', [0.41 0.955 0.08 0.035], 'String', call_items, ...
                  'Value', call_sel_idx, 'Callback', @thermoturb_coverage_filter_callback);
        info.call_popup  = call_popup;
        info.call_offset = call_offset;
    end

    % Second row: "Inner iterations" (1D-path logs only) and "Coolant
    % segment" (2D-path logs only) never both appear (see docstring), so
    % this row holds at most "Inner iterations" OR "Coolant segment"
    % alongside "Air segment" (which shows for either path) -- never all
    % three crowded together.
    if show_inner_dropdown
        uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
                  'Position', [0.01 0.91 0.1 0.03], 'String', 'Inner iterations:', ...
                  'HorizontalAlignment', 'left', 'BackgroundColor', fig.Color);
        inner_popup = uicontrol(fig, 'Style', 'popupmenu', 'Units', 'normalized', ...
                  'Position', [0.11 0.91 0.12 0.035], 'String', inner_items, ...
                  'Value', inner_sel_idx, 'Callback', @thermoturb_coverage_filter_callback);
        info.inner_popup  = inner_popup;
        info.inner_offset = inner_offset;
    end

    if has_air_seg
        uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
                  'Position', [0.25 0.91 0.1 0.03], 'String', 'Air segment:', ...
                  'HorizontalAlignment', 'left', 'BackgroundColor', fig.Color);
        air_seg_popup = uicontrol(fig, 'Style', 'popupmenu', 'Units', 'normalized', ...
                  'Position', [0.35 0.91 0.08 0.035], 'String', air_seg_items, ...
                  'Value', air_seg_sel_idx, 'Callback', @thermoturb_coverage_filter_callback);
        info.air_seg_popup = air_seg_popup;
        info.air_seg_rank  = air_seg_rank;
    end

    if show_coolant_dropdown
        uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
                  'Position', [0.46 0.91 0.13 0.03], 'String', 'Coolant segment:', ...
                  'HorizontalAlignment', 'left', 'BackgroundColor', fig.Color);
        coolant_seg_popup = uicontrol(fig, 'Style', 'popupmenu', 'Units', 'normalized', ...
                  'Position', [0.59 0.91 0.08 0.035], 'String', coolant_seg_items, ...
                  'Value', coolant_seg_sel_idx, 'Callback', @thermoturb_coverage_filter_callback);
        info.coolant_seg_popup = coolant_seg_popup;
        info.coolant_seg_rank  = coolant_seg_rank;
    end

    n_shown_init = sum(iter_mask_init & call_mask_init & inner_mask_init & air_seg_mask_init & coolant_seg_mask_init);
    info.count_text = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
              'Position', [0.73 0.955 0.25 0.03], ...
              'String', sprintf('%d / %d points shown', n_shown_init, n_pts), ...
              'HorizontalAlignment', 'left', 'BackgroundColor', fig.Color);

    setappdata(fig, 'coverage_info', info);
end

end

%% ---- Local helper: mask from an offset-based dropdown's selection ---- %%
% Shared by all three dropdowns (iter/call/inner), since all three use the
% same 'All' / 'Last (converged)' (offset 0) / '-1','-2',... scheme.
function mask = offset_mask_from_selection(offset_col, sel_str)
if strcmp(sel_str, 'All')
    mask = true(size(offset_col));
elseif strcmp(sel_str, 'Last (converged)')
    mask = offset_col == 0;
else
    mask = offset_col == str2double(sel_str);
end
end

%% ---- Local helper: an offset-based dropdown's item list for a row scope ---- %%
% 'All', 'Last (converged)' (offset 0), then whichever negative offsets
% actually occur among rows(scope_mask), most-recent-first (-1,-2,...) --
% e.g. once narrowed to one specific iteration (for the call-# dropdown)
% or one specific design run (for the inner dropdown), only as far back as
% that selection actually goes. 0 itself is never listed as a bare number
% since 'Last (converged)' already covers it.
function items = offset_items_for_scope(offset_col, scope_mask)
vals  = unique(offset_col(scope_mask));
vals  = sort(vals(vals ~= 0), 'descend');   % -1, -2, -3, ... (closest to last/converged first)
items = [{'All'; 'Last (converged)'}; cellstr(num2str(vals))];
end

%% ---- Local helper: dropdown index matching a desired offset value ---- %%
% Used only for the 'iter' function argument, which is still given as an
% absolute iteration number and converted to an offset by the caller.
% Returns [] if that offset isn't an available option.
function idx = offset_item_index(items, offset_val)
if offset_val == 0
    idx = find(strcmp(items, 'Last (converged)'), 1);
else
    idx = find(strcmp(items, sprintf('%d', offset_val)), 1);
end
end

%% ---- Local helper: mask from a forward-rank dropdown's selection ---- %%
% Shared by the Air segment / Coolant segment dropdowns, which filter on a
% PHYSICAL position (1 = first, 2 = second, ...) rather than a convergence
% loop, so they use plain forward numbering instead of the
% backward-from-converged 'Last (converged)'/-1/-2/... scheme above.
function mask = fwd_mask_from_selection(rank_col, sel_str)
if strcmp(sel_str, 'All')
    mask = true(size(rank_col));
else
    mask = rank_col == str2double(sel_str);
end
end

%% ---- Local helper: a forward-rank dropdown's item list for a row scope ---- %%
% 'All', then whichever ranks actually occur among rows(scope_mask), in
% ascending order (1,2,3,...) -- e.g. once narrowed to one mass-flow
% iteration, only as many segments/columns as that iteration's run
% actually had.
function items = fwd_items_for_scope(rank_col, scope_mask)
vals  = unique(rank_col(scope_mask));
vals  = sort(vals(vals > 0), 'ascend');
items = [{'All'}; cellstr(num2str(vals))];
end

%% ---- Local helper: update a popup's items, keeping the current ---- %%
%% ---- selection if it still exists, else falling back to default_idx ---- %%
function update_popup_items(popup, new_items, default_idx)
cur_items = popup.String;
cur_sel   = cur_items{popup.Value};
match     = find(strcmp(new_items, cur_sel), 1);
if isempty(match)
    match = default_idx;
end
set(popup, 'String', new_items, 'Value', match);
end

%% ---- Callback: mass-flow iteration / Q_pred_L call # / inner-iteration dropdowns ---- %%
% Shared by all three dropdowns (whichever one changed). The three filters
% combine (iter AND call AND inner), and each dropdown's available OPTIONS
% cascade from the coarser one(s): changing iter rebuilds which call
% offsets are offered (only those the selected iteration(s) actually
% reached), and changing either rebuilds which inner pass offsets are
% offered (only those the now-selected iter+call rows actually reached).
% The current selection in a dropdown being rebuilt is kept if it's still
% a valid option, otherwise it falls back to that dropdown's default
% ('All' for iter/call, 'Last (converged)' for inner).
function thermoturb_coverage_filter_callback(src, ~)
fig  = ancestor(src, 'figure');
info = getappdata(fig, 'coverage_info');

iter_items = info.iter_popup.String;
iter_sel   = iter_items{info.iter_popup.Value};
iter_mask  = offset_mask_from_selection(info.iter_offset, iter_sel);

if ~isempty(info.call_popup)
    new_call_items = offset_items_for_scope(info.call_offset, iter_mask);
    update_popup_items(info.call_popup, new_call_items, 1);   % 1 = 'All'
    call_items = info.call_popup.String;
    call_sel   = call_items{info.call_popup.Value};
    call_mask  = offset_mask_from_selection(info.call_offset, call_sel);
else
    call_mask = true(height(info.log_tbl), 1);
end

if ~isempty(info.inner_popup)
    new_inner_items = offset_items_for_scope(info.inner_offset, iter_mask & call_mask);
    update_popup_items(info.inner_popup, new_inner_items, 2);   % 2 = 'Last (converged)'
    inner_items = info.inner_popup.String;
    inner_sel   = inner_items{info.inner_popup.Value};
    inner_mask  = offset_mask_from_selection(info.inner_offset, inner_sel);
else
    inner_mask = true(height(info.log_tbl), 1);
end

% Air segment / Coolant segment are independent of each other (two
% different axes of the same k-by-l grid, neither nested in the other),
% so each is rescoped only by the coarser iter/call/inner selections, not
% by the other's current pick.
scope_for_seg = iter_mask & call_mask & inner_mask;

if ~isempty(info.air_seg_popup)
    new_air_seg_items = fwd_items_for_scope(info.air_seg_rank, scope_for_seg);
    update_popup_items(info.air_seg_popup, new_air_seg_items, 1);   % 1 = 'All'
    air_seg_items = info.air_seg_popup.String;
    air_seg_sel   = air_seg_items{info.air_seg_popup.Value};
    air_seg_mask  = fwd_mask_from_selection(info.air_seg_rank, air_seg_sel);
else
    air_seg_mask = true(height(info.log_tbl), 1);
end

if ~isempty(info.coolant_seg_popup)
    new_coolant_seg_items = fwd_items_for_scope(info.coolant_seg_rank, scope_for_seg);
    update_popup_items(info.coolant_seg_popup, new_coolant_seg_items, 1);   % 1 = 'All'
    coolant_seg_items = info.coolant_seg_popup.String;
    coolant_seg_sel   = coolant_seg_items{info.coolant_seg_popup.Value};
    coolant_seg_mask  = fwd_mask_from_selection(info.coolant_seg_rank, coolant_seg_sel);
else
    coolant_seg_mask = true(height(info.log_tbl), 1);
end

combined_mask = iter_mask & call_mask & inner_mask & air_seg_mask & coolant_seg_mask;

if ~isempty(info.count_text)
    set(info.count_text, 'String', sprintf('%d / %d points shown', sum(combined_mask), height(info.log_tbl)));
end

for p = 1:numel(info.panels)
    pan = info.panels(p);
    show_bulk = info.bulk_mask & combined_mask;
    set(pan.bulk_h, 'XData', pan.x_all(show_bulk), 'YData', pan.y_all(show_bulk), ...
                    'CData', info.color_idx(show_bulk));
    if strcmp(info.mode, '3d')
        set(pan.bulk_h, 'ZData', pan.z_all(show_bulk));
    end
    % Rescale the color span to whatever's now shown -- without this the
    % axes keeps whatever CLim the very first ('All') draw computed, so
    % filtering down to one iteration (or one call #) would otherwise
    % show only a sliver of the colormap (or, worse, look like a single
    % flat color if that selection's values happen to sit in a narrow
    % band of the original full-run range).
    if any(show_bulk)
        caxis(pan.ax, 'auto');
    end
    if ~isempty(pan.highlight_h)
        show_hl = info.is_highlight & combined_mask;
        set(pan.highlight_h, 'XData', pan.x_all(show_hl), 'YData', pan.y_all(show_hl));
        if strcmp(info.mode, '3d')
            set(pan.highlight_h, 'ZData', pan.z_all(show_hl));
        end
    end
end
end

%% ---- Local helper: narrowest-Pr-span Re range at one (Tb,Tw) corner ---- %%
function [re_lo, re_hi] = corner_re_range(tt, Tb_val, Tw_val, regime)
mask = tt.regime == regime & abs(tt.Tb - Tb_val) < 0.1 & abs(tt.Tw - Tw_val) < 0.1;
sub  = tt(mask, :);
Pr_vals = unique(sub.Pr);
n = numel(Pr_vals);
los = zeros(n,1); his = zeros(n,1); spans = zeros(n,1);
for j = 1:n
    rows   = abs(sub.Pr - Pr_vals(j)) < 1e-4;
    los(j) = min(sub.Re(rows));
    his(j) = max(sub.Re(rows));
    spans(j) = his(j) - los(j);
end
[~, sel] = min(spans);
re_lo = los(sel);
re_hi = his(sel);
end

%% ---- Local helper: one coverage level, hull ring projected onto Re vs. Tb or Tw (2D) ---- %%
% re_vals are the level's Re values at each hull corner, y_vals the
% matching Tb (or Tw) values; corners are in hull order, so closing the
% loop gives every hull edge projected onto the panel.
function draw_ring_lines(ax, re_vals, y_vals, col, ln, lw, name)
ring = [1:numel(re_vals) 1];
plot(ax, re_vals(ring), y_vals(ring), ln, 'Color', col, 'LineWidth', lw, 'HandleVisibility', 'off');
% Legend proxy.
plot(ax, nan, nan, ln, 'Color', col, 'LineWidth', lw, 'DisplayName', name);
end

%% ---- Local helper: 3D wireframe box through the 3 corner levels (3D) ---- %%
% corners rows are the (Tb,Tw) hull vertices in hull order (consecutive rows
% are neighbours). lo/mid/hi are each n x 1, the Re value of that level at
% each corner. Vertical edges run lo -> mid -> hi at each corner; each
% level also gets a ring around the hull. All three levels are
% wireframe-only -- no filled plane.
function draw_wireframe_box(ax, lo, mid, hi, corners, col, name)
n = size(corners, 1);
ring_order = [1:n 1];   % closed loop around the hull perimeter

plot3(ax, lo(ring_order),  corners(ring_order,1), corners(ring_order,2), '-',  'Color', col, 'LineWidth', 1.2, 'HandleVisibility', 'off');
plot3(ax, mid(ring_order), corners(ring_order,1), corners(ring_order,2), '--', 'Color', col, 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot3(ax, hi(ring_order),  corners(ring_order,1), corners(ring_order,2), '-',  'Color', col, 'LineWidth', 1.2, 'HandleVisibility', 'off');

for i = 1:n
    plot3(ax, [lo(i) mid(i) hi(i)], corners(i,1)*[1 1 1], corners(i,2)*[1 1 1], '-', ...
          'Color', col, 'LineWidth', 1.0, 'HandleVisibility', 'off');
end

plot3(ax, nan, nan, nan, '-', 'Color', col, 'LineWidth', 1.5, 'DisplayName', name);
end

%% ---- Local helper: dashed frame at a constant Re around the hull (3D) ---- %%
% Same look as the middle ("turbulent DNS data ends") ring of the
% wireframe box, but at a fixed Re. Used for the laminar/transition/
% turbulent regime boundaries (Re_lam_max, Re_tur_min).
function draw_re_frame(ax, re, corners, col, name)
n = size(corners, 1);
ring = [1:n 1];
plot3(ax, re*ones(1, n+1), corners(ring,1), corners(ring,2), '--', ...
      'Color', col, 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot3(ax, nan, nan, nan, '--', 'Color', col, 'LineWidth', 1.0, 'DisplayName', name);
end

%% ---- Local helper: 'on'/'off' for HandleVisibility from a logical ---- %%
function s = tern(tf)
if tf, s = 'on'; else, s = 'off'; end
end