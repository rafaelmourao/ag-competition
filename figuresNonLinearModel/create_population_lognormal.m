%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create_population_lognormal.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-03-18
%}

clear;
addpath('../classes')
rng(1);

nworkers = 90;

% Input model parameters
meanS = sqrt(25000^2 - 5100^2);
typeDistributionMean{1} = ...
    [1.9*1e-3, 1330, 4340, meanS]

typeDistributionLogCovariance{1} = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.98 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

costOfPublicFunds = 0;

% Calculation parameters
populationSize = 1e6
publicInsuranceMaximum = 30000

CalculationParametersEquilibrium.behavioralAgents = 0.1;
CalculationParametersEquilibrium.fudge            = 1e-4;
CalculationParametersEquilibrium.maxIterations    = 1e5;
CalculationParametersEquilibrium.tolerance        = 10;

CalculationParametersEquilibrium

CalculationParametersOptimum.maxIterations        = 10^3;
CalculationParametersOptimum.tolerance            = 0.01;

CalculationParametersOptimum

% Contracts

deductibleVector = [ publicInsuranceMaximum, 1500, 750, 500, 250, linspace(5200,0,20) ]
coinsuranceVector = [ 1, .1, .1, .1, .1, linspace(0.35,.1,20)  ]
oopMaxVector = [ publicInsuranceMaximum, 4500, 3750, 3500, 2750, linspace(6400,2500,20) ]

contracts = [ 1 6 14 18 21 23 24 25 ]

deductibleVector = deductibleVector(contracts);
coinsuranceVector = coinsuranceVector(contracts);
oopMaxVector = oopMaxVector(contracts);

% Calculate equilibrium and optimal prices.

if ~isempty(gcp('nocreate'))
     delete(gcp)
end

poolobj = parpool(nworkers)

for i = 1:length(typeDistributionMean)
    
Model(i) = healthcaralognormalmodel_nl( deductibleVector, ...
    coinsuranceVector, oopMaxVector, publicInsuranceMaximum, ...  
    typeDistributionMean{i}, typeDistributionLogCovariance{i});

for j = 1:length(deductibleVector)
    meanCoverage(j) = Model(i).meanCoverage(Model(i).contracts{j});
end

disp(meanCoverage)

tic
Population(i) = population(Model(i), populationSize, nworkers);
time=toc;
fprintf('Time to create sample: %f\n',time)

end

save Populations.mat

delete(poolobj)

