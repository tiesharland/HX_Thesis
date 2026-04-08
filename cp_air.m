function cp_air = cp_air(T)
cp_air = 1005 + 0.1*(T-288) - 0.0002*(T-288).^2;
end