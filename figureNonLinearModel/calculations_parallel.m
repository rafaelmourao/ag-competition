clear;
addpath('../classes');
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
populationSize = 10^6;

CalculationParametersEquilibrium.behavioralAgents = 0.01;
CalculationParametersEquilibrium.fudge            = 5 * 10^(-4);
CalculationParametersEquilibrium.maxIterations    = 10^5;
CalculationParametersEquilibrium.tolerance        = 10;

CalculationParametersOptimum.maxIterations        = 10^3;
CalculationParametersOptimum.tolerance            = 0.01;

% List of models
modelName{1}              = 'interval';
slopeVector{1}            = 0:0.04:1;
moralHazardLogVariance{1} = 0.28;

modelName{2}              = 'mandate';
slopeVector{2}            = 0.60:0.04:1;
moralHazardLogVariance{2} = 0.28;

modelName{3}              = 'interval_high_mh_variance';
slopeVector{3}            = 0:0.04:1;
moralHazardLogVariance{3} = 0.98;

modelName{4}              = 'mandate_high_mh_variance';
slopeVector{4}            = 0.60:0.04:1;
moralHazardLogVariance{4} = 0.98;

% Loop
nSimulations = length(modelName);

poolObject = parpool(nSimulations);

parfor i = 1 : nSimulations
    innerTypeDistributionLogCovariance = typeDistributionLogCovariance;

    innerTypeDistributionLogCovariance(2, 2) = moralHazardLogVariance{i};
    Model = healthcaralognormalmodel(slopeVector{i}, typeDistributionMean, innerTypeDistributionLogCovariance);
    
    Population = population(Model, populationSize);

    [pEquilibrium, DEquilibrium, ACEquilibrium, ComputationOutputEquilibrium] = ...
            Population.findequilibrium(CalculationParametersEquilibrium);
    WEquilibrium = Population.welfare(pEquilibrium, ...
                                          costOfPublicFunds);
    
    [pEfficient, WEfficient, ComputationOutputEfficient] = ...
            findefficient(Population, costOfPublicFunds, CalculationParametersOptimum);
    DEfficient = Population.demand(pEfficient);

    display(modelName{i});
    display(ComputationOutputEfficient);
    display(ComputationOutputEquilibrium);    
    parsave(modelName{i}, Model, pEquilibrium, DEquilibrium, ...
                ACEquilibrium, ComputationOutputEquilibrium, ...
            WEquilibrium, pEfficient, WEfficient, ...
            ComputationOutputEfficient, DEfficient, ...
            CalculationParametersEquilibrium, CalculationParametersOptimum, ...
            costOfPublicFunds, populationSize)


end;

delete(poolObject);