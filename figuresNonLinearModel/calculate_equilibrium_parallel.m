%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_equilibrium_parallel.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-03-18

In this script we calculate and store the results of the simulations with
the populations created with create_population_lognormal.m.

%}

clear;
addpath('../classes')
load('Populations.mat')
rng(1);

test(1).contracts = 1:Population(1).nContracts;
test(2).contracts = 1:Population(1).nContracts;
test(3).contracts = 1:Population(1).nContracts;

test(1).populations = 1;
test(2).populations = 2;
test(3).populations = 3;

[test.mandate] = deal(7); % index of the first mandate contract

[test.costOfPublicFunds] = deal(0);

CalculationParametersEquilibrium.behavioralAgents = 0.01;
CalculationParametersEquilibrium.fudge            = 1e-8;
CalculationParametersEquilibrium.maxIterations    = 1e4;
CalculationParametersEquilibrium.tolerance        = 1;
CalculationParametersEquilibrium.lineSearchErrorTolerance = 10;
CalculationParametersEquilibrium.lineSearchBeta   = .9;
display(CalculationParametersEquilibrium)

CalculationParametersOptimum.maxIterations        = 1e3;
CalculationParametersOptimum.tolerance            = 0.01;
CalculationParametersOptimum.knitro               = 'true';
CalculationParametersOptimum.knitroMultistartN    = 300;
display(CalculationParametersOptimum)

[test.CalculationParametersEquilibrium] = deal(CalculationParametersEquilibrium);
[test.CalculationParametersOptimum] = deal(CalculationParametersOptimum);

Population_old = Population;
iter = 0;

for i = 1:length(test)
    
    nPopulations = length(test(i).populations);
    
    for j = 1:nPopulations
        iter = iter + 1;
        CalculationParametersEquilibrium(iter) = test(i).CalculationParametersEquilibrium;
        CalculationParametersOptimum(iter) = test(i).CalculationParametersOptimum;
        costOfPublicFunds(iter) = test(i).costOfPublicFunds;
        Population(iter) = Population_old(test(i).populations(j));
        Population(iter).uMatrix = Population(iter).uMatrix(:,test(i).contracts);
        Population(iter).cMatrix = Population(iter).cMatrix(:,test(i).contracts);
        Population(iter).nContracts = length(test(i).contracts);
        test(i).Model(j) = Model(test(i).populations(j));
        test(i).Model(j).contracts = test(i).Model(j).contracts(:,test(i).contracts);
    end
    for j = (nPopulations+1):(2*nPopulations)
        iter = iter + 1;
        CalculationParametersEquilibrium(iter) = test(i).CalculationParametersEquilibrium;
        CalculationParametersOptimum(iter) = test(i).CalculationParametersOptimum;
        costOfPublicFunds(iter) = test(i).costOfPublicFunds;
        Population(iter) = Population(iter-nPopulations);
        Population(iter).uMatrix = Population(iter-nPopulations).uMatrix(:,test(i).mandate:end);
        Population(iter).cMatrix = Population(iter-nPopulations).cMatrix(:,test(i).mandate:end);
        Population(iter).nContracts = length(test(i).contracts)-(test(i).mandate-1);
        test(i).Model(j) = test(i).Model(j-nPopulations);
        test(i).Model(j).contracts = test(i).Model(j).contracts(:,test(i).mandate:end);
    end
    
    test(i).nPopulations = nPopulations;
    
end

nworkers = 2 * sum(cat(1,test.nPopulations));

if ~isempty(gcp('nocreate'))
    delete(gcp)
end

poolobj = parpool(nworkers);

parfor i = 1:nworkers
    
    [pEquilibrium{i}, DEquilibrium{i}, ACEquilibrium{i}, ComputationOutputEquilibrium{i}] = ...
        Population(i).findequilibrium(CalculationParametersEquilibrium(i))
    
    WEquilibrium{i} = Population(i).welfare(pEquilibrium{i}, costOfPublicFunds(i))
    
    [pEfficient{i}, WEfficient{i}, ComputationOutputEfficient{i}] = ...
        findefficient(Population(i), costOfPublicFunds(i), CalculationParametersOptimum(i), pEquilibrium{i})
    
    DEfficient{i} = Population(i).demand(pEfficient{i})
    
end

iter = 0;
for i = 1:length(test)
    for j = 1:2*test(i).nPopulations
        iter = iter + 1;
        test(i).pEquilibrium{j} = pEquilibrium{iter};
        test(i).DEquilibrium{j} = DEquilibrium{iter};
        test(i).ACEquilibrium{j} = ACEquilibrium{iter};
        test(i).WEquilibrium{j} = WEquilibrium{iter};
        test(i).ComputationOutputEquilibrium{j} = ComputationOutputEquilibrium{iter};
        test(i).pEfficient{j} = pEfficient{iter};
        test(i).WEfficient{j} = WEfficient{iter};
        test(i).ComputationOutputEfficient{j} = ComputationOutputEfficient{iter};
        test(i).DEfficient{j} = DEfficient{iter};
    end
end

delete(poolobj)
clear poolobj Population*
save('tests.mat')

print_lognormal_tables

