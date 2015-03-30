%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_equilibrium_parallel.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-03-18
%}

clear;
addpath('../../classes')
load('Populations.mat')
rng(1);

test(1).contracts = 1:Population(1).nContracts;
test(2).contracts = 1:Population(1).nContracts;
test(3).contracts = 1:Population(1).nContracts;
test(4).contracts = 1:Population(1).nContracts;

test(1).populations = 1:6;
test(2).populations = 7:12;
test(3).populations = 13:18;
test(4).populations = 19:21;

[test.mandate] = deal(7); % index of the first mandate contract

[test.costOfPublicFunds] = deal(0);

CalculationParametersEquilibrium.behavioralAgents = 0.01;
CalculationParametersEquilibrium.fudge            = 1e-6;
CalculationParametersEquilibrium.maxIterations    = 1e4;
CalculationParametersEquilibrium.tolerance        = 1;
CalculationParametersEquilibrium.lineSearchErrorTolerance = Inf;
display(CalculationParametersEquilibrium)

CalculationParametersOptimum.maxIterations        = 1e3;
CalculationParametersOptimum.tolerance            = 0.01;
CalculationParametersOptimum.knitro               = 'true';
display(CalculationParametersOptimum)

[test.CalculationParametersEquilibrium] = deal(CalculationParametersEquilibrium);
[test.CalculationParametersOptimum] = deal(CalculationParametersOptimum);

Population_old = Population;
iter = 0;

for i = 1:length(test)

    nPopulations = length(test(i).populations);
    
for j = 1:nPopulations
    iter = iter + 1;
    CalculationParametersEquilibrium(i) = test(i).CalculationParametersEquilibrium;
    CalculationParametersOptimum(i) = test(i).CalculationParametersOptimum;
    costOfPublicFunds(i) = test(i).costOfPublicFunds;    
    Population(iter) = Population_old(test(i).populations(j));
    Population(iter).uMatrix = Population(iter).uMatrix(:,test(i).contracts);
    Population(iter).cMatrix = Population(iter).cMatrix(:,test(i).contracts);
    Population(iter).nContracts = length(test(i).contracts);
    test(i).Model(j) = Model(test(i).populations(j));
    test(i).Model(j).contracts = test(i).Model(j).contracts(:,test(i).contracts);
end
for j = (nPopulations+1):(2*nPopulations)
    iter = iter + 1;
    CalculationParametersEquilibrium(i) = test(i).CalculationParametersEquilibrium;
    CalculationParametersOptimum(i) = test(i).CalculationParametersOptimum;
    costOfPublicFunds(i) = test(i).costOfPublicFunds;    
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

poolobj = parpool(nworkers)

parfor i = 1:nworkers
    
[pEquilibrium{i}, DEquilibrium{i}, ACEquilibrium{i}, ComputationOutputEquilibrium{i}] = ...
         Population(i).findequilibrium(CalculationParametersEquilibrium(i))
WEquilibrium{i} = Population(i).welfare(pEquilibrium{i}, costOfPublicFunds(i))

[pEfficient{i}, WEfficient{i}, ComputationOutputEfficient{i}] = ...
            findefficient(Population(i), costOfPublicFunds(i), CalculationParametersOptimum(i), pEquilibrium{i})
        
DEfficient{i} = Population(i).demand(pEfficient{i})

end

for i = 1:length(test)
    for j = 1:test(i).nPopulations
    z = (2*nPopulations)*(i-1) + j;
    test(i).pEquilibrium{j} = pEquilibrium{z};
    test(i).DEquilibrium{j} = DEquilibrium{z};
    test(i).ACEquilibrium{j} = ACEquilibrium{z};
    test(i).WEquilibrium{j} = WEquilibrium{z};
    test(i).ComputationOutputEquilibrium{j} = ComputationOutputEquilibrium{z};
    test(i).pEfficient{j} = pEfficient{z};
    test(i).WEfficient{j} = WEfficient{z};
    test(i).ComputationOutputEfficient{j} = ComputationOutputEfficient{z};
    test(i).DEfficient{j} = DEfficient{z};
    end
end

delete(poolobj)
clear poolobj Population*
save('tests.mat')

print_lognormal_tables

