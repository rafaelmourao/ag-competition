%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_caralognormal_model_nl_equilibrium %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-03-02

This script is used to test the equilibrium calculations for the 
healthcaralognormalmodel_nl class

%}

clear;
addpath('../classes')
rng(1);

% Input model parameters
meanS = sqrt(25000^2 - 5100^2);
typeDistributionMean = ...
    [1.5*1e-3, 1330, 4340, meanS]; % Original A was 1.9*10^-3
typeDistributionLogCovariance = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25]; % ???

% Calculation parameters
populationSize = 1e6;
publicInsuranceMaximum = 30000;

CalculationParametersEquilibrium.behavioralAgents = 0.1;
CalculationParametersEquilibrium.fudge            = 1e-5;
CalculationParametersEquilibrium.maxIterations    = 1e5;
CalculationParametersEquilibrium.tolerance        = 10;



% Contracts

deductibleVector = [ publicInsuranceMaximum, 1500, 750, 500, 250, 0 ];
coinsuranceVector = [ 1, .1, .1, .1, .1, .1 ];
oopMaxVector = [ publicInsuranceMaximum, 4500, 3750, 3500, 2750, 2500 ];


% Calculate equilibrium and optimal prices.
Model = healthcaralognormalmodel_nl( deductibleVector, ...
    coinsuranceVector, oopMaxVector, publicInsuranceMaximum, ...  
    typeDistributionMean, typeDistributionLogCovariance);

if ~isempty(gcp)
     delete(gcp)
end

parpool(32)
tic
Population = population(Model, populationSize, 32);
time=toc;
delete(gcp)

fprintf('Time to create sample: %f\n',time)

[pEquilibrium, DEquilibrium, ACEquilibrium, ComputationOutputEquilibrium] = ...
         Population.findequilibrium(CalculationParametersEquilibrium)

clear Population
save equilibrium_results.mat