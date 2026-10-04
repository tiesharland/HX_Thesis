function out = thermoturb_query_log(mode, entry)
%THERMOTURB_QUERY_LOG  Persistent log of every thermoturb_cached query,
%so you can check afterward how much of a model run's DNS lookups fell
%inside the table's coverage vs. were extrapolated/clamped.
%
%   thermoturb_query_log('log', entry)   -- append one query. Called by
%                                            thermoturb_cached itself;
%                                            you shouldn't need to call
%                                            this directly.
%   log = thermoturb_query_log('get')    -- return everything logged so
%                                            far, as a table with columns
%                                            Re, Pr, Tb, Tw, regime,
%                                            Re_extrap, Pr_extrap,
%                                            Tb_extrap, Tw_extrap,
%                                            any_extrap.
%   thermoturb_query_log('reset')        -- clear the log. Call this
%                                            before a run if you only
%                                            want that run's queries
%                                            (otherwise it accumulates
%                                            across every thermoturb_cached
%                                            call made this session).
%
% entry (for 'log'): a 1x1 struct with scalar fields Re, Pr, Tb, Tw,
% regime (0 = laminar, 1 = turbulent, 2 = transition blend of both), and
% the four logical *_extrap fields -- true means that dimension fell
% outside the DNS table's coverage for whichever regime(s) were actually
% evaluated for this query.
%
% Storage is preallocated and doubled on growth (not a table that grows
% row-by-row), since thermoturb_cached may call this on every segment of
% every model run.

persistent Re_log Pr_log Tb_log Tw_log regime_log ...
    Re_extrap_log Pr_extrap_log Tb_extrap_log Tw_extrap_log n cap

if isempty(n)
    n   = 0;
    cap = 1024;
    Re_log        = zeros(cap, 1);
    Pr_log        = zeros(cap, 1);
    Tb_log        = zeros(cap, 1);
    Tw_log        = zeros(cap, 1);
    regime_log    = zeros(cap, 1);
    Re_extrap_log = false(cap, 1);
    Pr_extrap_log = false(cap, 1);
    Tb_extrap_log = false(cap, 1);
    Tw_extrap_log = false(cap, 1);
end

switch mode
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
        any_extrap = Re_extrap | Pr_extrap | Tb_extrap | Tw_extrap;
        out = table(Re, Pr, Tb, Tw, regime, Re_extrap, Pr_extrap, Tb_extrap, Tw_extrap, any_extrap);

    case 'reset'
        n   = 0;
        out = [];

    otherwise
        error('thermoturb_query_log:UnknownMode', 'Unknown mode "%s".', mode);
end

end