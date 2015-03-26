%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_equilibrium_parallel.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-03-18
%}

clear;
addpath('../classes')
load('Populations.mat')
rng(1);


costOfPublicFunds = 0;
CalculationParametersEquilibrium.behavioralAgents = 0.01;
CalculationParametersEquilibrium.fudge            = 5e-5;
CalculationParametersEquilibrium.maxIterations    = 1e6;
CalculationParametersEquilibrium.tolerance        = 1;
display(CalculationParametersEquilibrium)
CalculationParametersOptimum.maxIterations        = 1e3;
CalculationParametersOptimum.tolerance            = 0.01;
display(CalculationParametersOptimum)

nPopulations = length(Population);
nworkers = 2*length(Population) * length(test);

Population_old = Population;
iter = 0;

for i = 1:length(test) 
    
test(i).costOfPublicFunds = costOfPublicFunds;
test(i).CalculationParametersEquilibrium = CalculationParametersEquilibrium;
test(i).CalculationParametersOptimum = CalculationParametersOptimum;
for j = 1:nPopulations
    iter = iter + 1;
    Population(iter) = Population_old(j);
    Population(iter).uMatrix = Population(iter).uMatrix(:,test(i).contracts);
    Population(iter).cMatrix = Population(iter).cMatrix(:,test(i).contracts);
    Population(iter).nContracts = length(test(i).contracts);
    test(i).Model(j) = Model(j);
    test(i).Model(j).contracts = test(i).Model(j).contracts(:,test(i).contracts);
end
for j = (nPopulations+1):(2*nPopulations)
    iter = iter + 1;
    Population(iter) = Population(iter-nPopulations);
    Population(iter).uMatrix = Population(iter-nPopulations).uMatrix(:,2:end);
    Population(iter).cMatrix = Population(iter-nPopulations).cMatrix(:,2:end);
    Population(iter).nContracts = length(test(i).contracts)-1;
    test(i).Model(j) = test(i).Model(j-nPopulations);
    test(i).Model(j).contracts = test(i).Model(j).contracts(:,2:end);
end

end

if ~isempty(gcp('nocreate'))
     delete(gcp)
end
poolobj = parpool(nworkers)

parfor i = 1:nworkers
    
[pEquilibrium{i}, DEquilibrium{i}, ACEquilibrium{i}, ComputationOutputEquilibrium{i}] = ...
         Population(i).findequilibrium(CalculationParametersEquilibrium)
WEquilibrium{i} = Population(i).welfare(pEquilibrium{i}, costOfPublicFunds)

[pEfficient{i}, WEfficient{i}, ComputationOutputEfficient{i}] = ...
            findefficient(Population(i), costOfPublicFunds, CalculationParametersOptimum)
        
DEfficient{i} = Population(i).demand(pEfficient{i})

end

for i = 1:length(test)
    for j = 1:(2*nPopulations)
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

clear Population*
save('tests.mat')

print_lognormal_tables

