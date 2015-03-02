%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_caralognormal_model_nl %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-03-02

This script simulates a CARA lognormal health nonlinear insurance model,
similar to that used in the Azevedo and Gottlieb 2015 paper. The simulation is run
with low precision, and not used in the paper. The purpose of this file is
to test the implementation of the class test_caralognormal_model_nl as well
as alternative functions
%}


clear;
addpath('../classes')
rng(1);

% Input model parameters
meanS = sqrt(25000^2 - 5100^2);
typeDistributionMean = ...
    [1*10^(-5), 1330, 4340, meanS]; % Original A was 1.9*10^-3
typeDistributionLogCovariance = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25]; % ???

costOfPublicFunds = 0;

% Calculation parameters
populationSize = 10^3;

CalculationParametersEquilibrium.behavioralAgents = 0.1;
CalculationParametersEquilibrium.fudge            = 1e-4;
CalculationParametersEquilibrium.maxIterations    = 5e4;
CalculationParametersEquilibrium.tolerance        = 0.01;

CalculationParametersOptimum.maxIterations        = 1e3;
CalculationParametersOptimum.tolerance            = 0.001;

% Contracts
deductibleVector = [ 1500, 750, 500, 250, 0 ];
coinsuranceVector = [ .1, .1, .1, .1, .1 ];
oopMaxVector = [ 4500, 3750, 3500, 2750, 2500 ];

% Calculate equilibrium and optimal prices.
Model = healthcaralognormalmodel_nl(deductibleVector, ...
    coinsuranceVector, oopMaxVector, typeDistributionMean, typeDistributionLogCovariance);

n = 1000;
u = zeros(1,n);
u_alt = zeros(1,n);
type=cell(1,n);
for i = 1:n
    type{i} = Model.typeDistribution();
end
for j = 1:3
    tic
    for i = 1:n
        u(i) = Model.uFunction(Model.contracts{1},type{i});
    end
    time(j)=toc;
end
disp('First test: Integrals by interval, utilities by interval')
fprintf('Sum of utilities: %f, Calculation times: %f %f %f seconds\n',sum(u),time)
for j = 1:3
    tic
    for i = 1:n
        u_alt(i) = Model.uFunction_alt(Model.contracts{1},type{i});
    end
    time(j)=toc;
end
disp('Second test: Integrals by whole positive region, utilities by interval')
fprintf('Sum of utilities: %f, Calculation time: Calculation times: %f %f %f seconds\n',sum(u),time)
for j = 1:3
    tic
    for i = 1:n
        u_alt2(i) = Model.uFunction_alt2(Model.contracts{1},type{i});
    end
    time(j)=toc;
end
disp('Third test: Integrals by interval, utilities by maximum over three cases')
fprintf('Sum of utilities: %f, Calculation time: Calculation times: %f %f %f seconds\n',sum(u),time)
for j = 1:3
    tic
    for i = 1:n
        u_alt3(i) = Model.uFunction_alt3(Model.contracts{1},type{i});
    end
    time(j)=toc;
end
disp('Fourth test: Integrals by whole positive region, utilities by maximum over three cases')
fprintf('Sum of utilities: %f, Calculation time: Calculation times: %f %f %f seconds\n',sum(u),time)
tic
Population = population(Model, populationSize);
time=toc;
fprintf('Time to create sample: %f\n',time)
[pEquilibrium, DEquilibrium, ACEquilibrium, ComputationOutputEquilibrium] = ...
         Population.findequilibrium(CalculationParametersEquilibrium)
% WEquilibrium = Population.welfare(pEquilibrium, ...
%                                   costOfPublicFunds);
%
% [pEfficient, WEfficient, ComputationOutputEfficient] = ...
%         findefficient(Population, costOfPublicFunds, CalculationParametersOptimum);
% DEfficient = Population.demand(pEfficient);