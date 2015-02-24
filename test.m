clear;
close all;

slopeVector = [0.00 0.90];
% slopeVector = 0:0.5:1;

meanS = sqrt(25000^2 - 5100^2);
typeDistributionMean = ...
    [1*10^(-5), 1330, 4340, meanS]; % Original A was 1.9*10^-3
typeDistributionLogCovariance = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
     -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
     -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
      0     0     0    0.01]; % ???
  
costOfPublicFunds = 0;

% populationSize = 1000;  
% CalculationParametersEquilibrium.behavioralAgents = 0.1;
% CalculationParametersEquilibrium.fudge            = 0.01;
% CalculationParametersEquilibrium.maxIterations    = 2000;
% CalculationParametersEquilibrium.tolerance        = 50;
% 
% CalculationParametersOptimum.maxIterations = 2000;
% CalculationParametersOptimum.tolerance     = 0.1;
%   
% Model = healthcaralognormalmodel(slopeVector, typeDistributionMean, typeDistributionLogCovariance);
% Population = population(Model, populationSize);

% Test basic functions
% p = linspace(0, 8000, Model.nContracts);
% [D, TC, CS, ~] = Population.demand(p);
% [p, D, AC, ComputationOutput] = Population.findequilibrium(CalculationParametersEquilibrium);

% Test computational functions
% [p, W, ComputationOutput] = findefficient(Population, costOfPublicFunds, CalculationParametersOptimum)
% Population.demand(p)

% Test graphing functions
% Population.graphEFC();

% % Test suggested calculation parameters
% [populationSize, CalculationParametersEquilibrium, CalculationParametersOptimum] = ...
%     Model.suggestComputationParameters(0.2);

% Test computational engine
clear slopeVector Model Population
slopeVector = 0:0.2:1;
Model = healthcaralognormalmodel(slopeVector, typeDistributionMean, typeDistributionLogCovariance);
% CalculationParametersEquilibrium.tolerance = CalculationParametersEquilibrium.tolerance / 10;
populationSize = 10^4;
CalculationParametersEquilibrium.behavioralAgents = 0.1;
CalculationParametersEquilibrium.tolerance        = 50;
CalculationParametersEquilibrium.maxIterations    = 10^3;
% CalculationParametersEquilibrium.fudge            = 10^(-2);

Population = population(Model, populationSize);
populationSize
CalculationParametersEquilibrium

[p, D, AC, ComputationOutput] = Population.findequilibriumalt(CalculationParametersEquilibrium);
ComputationOutput
p