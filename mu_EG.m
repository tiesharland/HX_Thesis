function mu = mu_EG(T)
    % Viscosity of ethylene glycol (Pa·s)
    mu = (2.2705e07*exp(-0.04675*T) + 1.158)*1e-03;
end