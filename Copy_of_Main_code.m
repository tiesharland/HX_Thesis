clear; clc;

%% Design variables for optimization %%
vars = [ ...  
    ... optimizableVariable('fr_dia',[0.4 0.6]), ...
    optimizableVariable('air_side',[2 12]), ...
    optimizableVariable('cool_side',[3.5 4.5]) ...
    ... optimizableVariable('AR_diff',[3 4]), ...
    ... optimizableVariable('AR_noz',[0.3 0.4]), ...
    ... optimizableVariable('N_segments', [5 100]) ...  
    ];

    % optimizableVariable('cool_side',[2 5]), ...


p = gcp('nocreate');   

if isempty(p)
    parpool(4);      
end

%% Bayesian Optimization %%
results = bayesopt(@Copy_of_tmsObjective, vars, 'MaxObjectiveEvaluations', 200,'IsObjectiveDeterministic', true,'AcquisitionFunctionName','expected-improvement-plus','Verbose',1,'UseParallel',true);

%% Get optimal configuration %%
bestConfig = results.XAtMinEstimatedObjective;

% fprintf(' AR_dif = %.3f, frontal dia = %.3f, air side fin = %f', bestConfig.AR_diff, bestConfig.fr_dia, bestConfig.air_side);
fprintf(' cool side fin = %.3f, air side fin = %f', bestConfig.cool_side, bestConfig.air_side);

