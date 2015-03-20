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

test(1).contracts = [ 1 2 3 4 5 25 ];
test(2).contracts = [ 1 6 14 18 21 23 24 25 ];

costOfPublicFunds = 0;
CalculationParametersEquilibrium.behavioralAgents = 0.1;
CalculationParametersEquilibrium.fudge            = 1e-4;
CalculationParametersEquilibrium.maxIterations    = 1e5;
CalculationParametersEquilibrium.tolerance        = 10;
CalculationParametersEquilibrium
CalculationParametersOptimum.maxIterations        = 1e3;
CalculationParametersOptimum.tolerance            = 0.01;
CalculationParametersOptimum

nworkers = 8 * length(test);

for i = 1:length(test) 

test(i).costOfPublicFunds = costOfPublicFunds;
test(i).CalculationParametersEquilibrium = CalculationParametersEquilibrium;
test(i).CalculationParametersOptimum = CalculationParametersOptimum;

for j = 1:4
    test(i).Population(j) = Population(j);
    test(i).Population(j).uMatrix = test(i).Population(j).uMatrix(:,test(i).contracts);
    test(i).Population(j).cMatrix = test(i).Population(j).cMatrix(:,test(i).contracts);
    test(i).Population(j).nContracts = length(test(i).contracts);
    test(i).Model(j) = Model(j);
    test(i).Model(j).contracts = test(i).Model(j).contracts(:,test(i).contracts);
end
for j = 5:8
    test(i).Population(j) = test(i).Population(j-4);
    test(i).Population(j).uMatrix = test(i).Population(j-4).uMatrix(:,2:end);
    test(i).Population(j).cMatrix = test(i).Population(j-4).cMatrix(:,2:end);
    test(i).Population(j).nContracts = length(test(i).contracts)-1;
    test(i).Model(j) = test(i).Model(j-4);
    test(i).Model(j).contracts = test(i).Model(j).contracts(:,2:end);
end

end

if ~isempty(gcp('nocreate'))
     delete(gcp)
end
poolobj = parpool(nworkers)

for i = 1:length(test)
    for j = 1:8
    z = 8*(i-1) + j;
    Population(z) = test(i).Population(j);
    end
end

parfor i = 1:nworkers
    
[pEquilibrium{i}, DEquilibrium{i}, ACEquilibrium{i}, ComputationOutputEquilibrium{i}] = ...
         Population(i).findequilibrium(CalculationParametersEquilibrium)


[pEfficient{i}, WEfficient{i}, ComputationOutputEfficient{i}] = ...
            findefficient(Population(i), costOfPublicFunds, CalculationParametersOptimum)
        
DEfficient{i} = Population(i).demand(pEfficient{i})

end

for i = 1:length(test)
    for j = 1:8
    z = 8*(i-1) + j;
    test(i).pEquilibrium{j} = pEquilibrium{z};
    test(i).DEquilibrium{j} = DEquilibrium{z};
    test(i).ACEquilibrium{j} = ACEquilibrium{z};
    test(i).ComputationOutputEquilibrium{j} = ComputationOutputEquilibrium{z};
    test(i).pEfficient{j} = pEfficient{z};
    test(i).WEfficient{j} = WEfficient{z};
    test(i).ComputationOutputEfficient{j} = ComputationOutputEfficient{z};
    test(i).DEfficient{j} = DEfficient{z};
    end
end

clear Population

save tests.mat

delete(poolobj)

