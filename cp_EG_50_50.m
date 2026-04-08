function cp = cp_EG_50_50(T)

    a3 = 4.23387495e-8;
    a2 = -4.15801375e-6;
    a1 = 3.06561504;
    a0 = 3.41161599e3;

    % Convert Kelvin to Celsius
    Tc = T - 273.15;
    
    cp = a3 .* Tc.^3 + a2 .* Tc.^2 + a1 .* Tc + a0; %digitzed from https://corecheminc.com/ethylene-glycol-water-mixture-properties/
end

