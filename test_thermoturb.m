%% Standalone test of thermoturb_query
% Compare against the known-good values from manual browser test:
% Re_bulk=3532, Pr=0.71, Tb=264, Tw=358
% Expected (from coefficients.dat, last converged row):
%   Rebulk  ~ 3532.0
%   Cf      ~ 4.7545e-3
%   Nusselt ~ 5.4004

clear
clc

Re_bulk = 3532;
Pr      = 0.71;
Tb      = 264;
Tw      = 358;

p = gcp('nocreate');   

if isempty(p)
    parpool(4);      
end

tic
[Re_out, Nu_out, Cf_out] = thermoturb_query(Re_bulk, Pr, Tb, Tw);
elapsed = toc;

fprintf('--- thermoturb_query result ---\n');
fprintf('Re_out  = %.4f   (expected ~3532)\n', Re_out);
fprintf('Nu_out  = %.4f   (expected ~5.4004)\n', Nu_out);
fprintf('Cf_out  = %.6e (expected ~4.7545e-3)\n', Cf_out);
fprintf('Elapsed time: %.3f s\n', elapsed);

%% Second call — test whether cookie/session handling works repeatedly
tic
[Re_out2, Nu_out2, Cf_out2] = thermoturb_query(1500, 0.71, 280, 350);
elapsed2 = toc;

fprintf('\n--- second call (different inputs) ---\n');
fprintf('Re_out2 = %.4f\n', Re_out2);
fprintf('Nu_out2 = %.4f\n', Nu_out2);
fprintf('Cf_out2 = %.6e\n', Cf_out2);
fprintf('Elapsed time: %.3f s\n', elapsed2);