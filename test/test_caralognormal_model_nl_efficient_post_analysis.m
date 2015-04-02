%% Initialize
addpath('../classes');
clear;
close all;
load test_caralognormal_model_nl_efficient.mat;

WEquilibrium = Population.welfare(pEquilibrium, 0);

DEfficient = Population.demand(pEfficient);

%% Plot graph
f = @(pvar) WEquilibrium - Population.welfare(pvar, 0);
n = 100;

alphagrid = linspace(-0.2, 1, n);
for i = 1:n
    p = alphagrid(i) .* pEfficient + (1-alphagrid(i)) .* pEquilibrium;
    fgrid(i) = f(p);
end;

plot(alphagrid, fgrid);


%% Optimize again
% Population = population(Model, 5, -1);
[pEfficient2, WEfficient2, ComputationOutputEfficient2] = ...
    Population.findefficient(costOfPublicFunds, CalculationParametersOptimum, pEquilibrium);