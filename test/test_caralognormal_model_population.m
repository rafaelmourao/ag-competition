%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_caralognormal_model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-03-4

This script tests performance for demand calculation for the
population class.
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

CalculationParametersEquilibrium.behavioralAgents = 0.2;
CalculationParametersEquilibrium.fudge            = 5 * 10^(-4);
CalculationParametersEquilibrium.maxIterations    = 10^3;
CalculationParametersEquilibrium.tolerance        = 1e-6;

results=struct;

for populationSize = [1e4, 1e5, 1e6]
    for ncontracts = [5 10 20 30 50 100 200]
        
        modelName              = 'interval';
        slopeVector            = linspace(0,1,ncontracts);
        
        disp(populationSize)
        disp(ncontracts)
        
        % Calculate equilibrium and optimal prices.
        Model = healthcaralognormalmodel(slopeVector, typeDistributionMean, typeDistributionLogCovariance);
        
        tic
        Population = population(Model, populationSize);        
        tempo = toc;
        fprintf('Population generation took %f seconds\n',tempo)
        
        disp(populationSize)
        disp(ncontracts)
        
        tic
        [pEquilibrium, DEquilibrium, ACEquilibrium, ComputationOutputEquilibrium] = ...
            Population.findequilibrium(CalculationParametersEquilibrium);
        toc
        
        results.(['p' num2str(populationSize)]).(['c' num2str(ncontracts)]).('t1')= ComputationOutputEquilibrium;
        
        tic
        [pEquilibrium, DEquilibrium, ACEquilibrium, ComputationOutputEquilibrium] = ...
            Population.findequilibrium_2(CalculationParametersEquilibrium);
        toc
        
        results.(['p' num2str(populationSize)]).(['c' num2str(ncontracts)]).('t2')= ComputationOutputEquilibrium;
        
        tic
        [pEquilibrium, DEquilibrium, ACEquilibrium, ComputationOutputEquilibrium] = ...
            Population.findequilibrium_3(CalculationParametersEquilibrium);
        toc
        
        results.(['p' num2str(populationSize)]).(['c' num2str(ncontracts)]).('t3')= ComputationOutputEquilibrium;
        
                tic
        [pEquilibrium, DEquilibrium, ACEquilibrium, ComputationOutputEquilibrium] = ...
            Population.findequilibrium_4(CalculationParametersEquilibrium);
        toc
        
        results.(['p' num2str(populationSize)]).(['c' num2str(ncontracts)]).('t4')= ComputationOutputEquilibrium;
        
                tic
        [pEquilibrium, DEquilibrium, ACEquilibrium, ComputationOutputEquilibrium] = ...
            Population.findequilibrium_5(CalculationParametersEquilibrium);
        toc
        
        results.(['p' num2str(populationSize)]).(['c' num2str(ncontracts)]).('t5')= ComputationOutputEquilibrium;
        
                        tic
        [pEquilibrium, DEquilibrium, ACEquilibrium, ComputationOutputEquilibrium] = ...
            Population.findequilibrium_6(CalculationParametersEquilibrium);
        toc
        
        results.(['p' num2str(populationSize)]).(['c' num2str(ncontracts)]).('t6')= ComputationOutputEquilibrium;
       
        save test_population_results.mat results
    end
end