function [value, isterminal, direction] = choke_event(~, M)
        value = 1 - M;   % stop when M → 1
        isterminal = 1;  % terminate integration
        direction = -1;  % only when approaching from subsonic side
end