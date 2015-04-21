%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_caralognormal_model_nl_efficient %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-03-02

This script is used to test the efficiency calculations for the 
healthcaralognormalmodel_nl class

%}

clear;
addpath('../classes')
rng(1);

% Input model parameters
meanS = sqrt(25000^2 - 5100^2);
typeDistributionMean = ...
    [1.9*1e-3, 1330, 4340, meanS]; % Original A was 1.9*10^-3
typeDistributionLogCovariance = ...
     [  0.25 -0.00 -0.00 0     ;
     -0.00   0.01 -0.00 0     ;
     -0.00 -0.00   0.01 0     ;
     0        0        0      0.01];


costOfPublicFunds = 0;

% Calculation parameters
populationSize = 1e5;
publicInsuranceMaximum = 30000;

CalculationParametersOptimum.maxIterations        = 10^3;
CalculationParametersOptimum.tolerance            = 0.01;

CalculationParametersEquilibrium.behavioralAgents = 0.1;
CalculationParametersEquilibrium.fudge            = 1e-4;
CalculationParametersEquilibrium.maxIterations    = 1e5;
CalculationParametersEquilibrium.tolerance        = 10;

nWorkers = 90;

% Contracts
deductibleVector = [ publicInsuranceMaximum, 1500, 750, 500, 250, linspace(5200,0,20) ];
coinsuranceVector = [ 1, .1, .1, .1, .1, linspace(0.35,.1,20)  ];
oopMaxVector = [ publicInsuranceMaximum, 4500, 3750, 3500, 2750, linspace(6400,2500,20) ];

contracts = [ 1 6 14 18 21 23 24 25 ];

deductibleVector = deductibleVector(contracts);
coinsuranceVector = coinsuranceVector(contracts);
oopMaxVector = oopMaxVector(contracts);

% Set up model
Model = healthcaralognormalmodel_nl( deductibleVector, ...
    coinsuranceVector, oopMaxVector, publicInsuranceMaximum, ...  
    typeDistributionMean, typeDistributionLogCovariance);


% Calculate population
if ~isempty(gcp('nocreate'))
     delete(gcp)
end

poolobj = parpool(nWorkers)

Population = population(Model, populationSize, nWorkers);

delete(poolobj)

% Calculate optimal and efficient prices.
[pEfficient, WEfficient, ComputationOutputEfficient] = Population.findefficient(costOfPublicFunds, CalculationParametersOptimum);
[pEquilibrium, DEquilibrium, ACEquilibrium, ComputationOutputEquilibrium] = ...
         Population.findequilibrium(CalculationParametersEquilibrium);

save test_caralognormal_model_nl_efficient.mat;