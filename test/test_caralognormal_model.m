%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_caralognormal_model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-02-24

This script simulates a CARA lognormal health insurance model, similar to
that used in the Azevedo and Gottlieb 2015 paper. The simulation is run
with low precision, and not used in the paper. The purpose of this file is
to test the classes in the agSim directory.
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
populationSize = 1e3;

CalculationParametersEquilibrium.behavioralAgents = 0.2;
CalculationParametersEquilibrium.fudge            = 1e-6;
CalculationParametersEquilibrium.maxIterations    = 1e4;
CalculationParametersEquilibrium.tolerance        = 1;

CalculationParametersOptimum.maxIterations        = 1e3;
CalculationParametersOptimum.tolerance            = 0.01;

% List of models
modelName              = 'interval';
slopeVector            = 0:0.2:1;

% Calculate equilibrium prices.
Model = healthcaralognormalmodel(slopeVector, typeDistributionMean, typeDistributionLogCovariance);
Population = population(Model, populationSize);

[pEquilibrium, DEquilibrium, ACEquilibrium, ComputationOutputEquilibrium] = ...
        Population.findequilibrium(CalculationParametersEquilibrium);
WEquilibrium = Population.welfare(pEquilibrium, ...
                                  costOfPublicFunds);

display(ComputationOutputEquilibrium);
display(pEquilibrium);
display(DEquilibrium);
display(WEquilibrium);

% Efficient prices with the original algorithm.                              
[pEfficient, WEfficient, ComputationOutputEfficient] = ...
        findefficient(Population, costOfPublicFunds, CalculationParametersOptimum, pEquilibrium);
DEfficient = Population.demand(pEfficient);

display(ComputationOutputEfficient);
display(pEfficient);
display(DEfficient);
display(WEfficient);

% Efficient prices with knitro.
CalculationParametersOptimum.knitro = 'true';
[pEfficient, WEfficient, ComputationOutputEfficient] = ...
        findefficient(Population, costOfPublicFunds, CalculationParametersOptimum, pEquilibrium);
DEfficient = Population.demand(pEfficient);

display(ComputationOutputEfficient);
display(pEfficient);
display(DEfficient);
display(WEfficient);