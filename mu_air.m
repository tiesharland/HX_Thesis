function mu_air = mu_air (T)

% Sutherland constants for air
  C = 1.458e-6;   % [kg/(m·s·√K)]
  S = 110.4;      % [K]

% Sutherland's formula
   mu_air = C .* (T.^1.5) ./ (T + S);

%mu_air = 4E-29*h^5 - 6E-24*h^4 + 1E-19*h^3 + 1E-14*h^2 - 4E-10*h + 2E-05; 
end