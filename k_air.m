function k_air = k_air(T)
% Thermal conductivity of air (W/m·K)
k_air = (1.5207e-11)*T.^3 - (4.8574e-08)*T.^2 + (1.0184e-04)*T - 3.9333e-04;
%k_air = 0.0177 + (4.03e-05)*T + (-1e-08)*T^2; %NASA correlation
end
