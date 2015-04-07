clear;
addpath('../classes');
load('Populations.mat','Population')
load('tests.mat','test')

for i = 1:3

populationSize           = Population(i).size;
nPopulationDemandProfile = min(7000, populationSize);

Interval.Model = test(i).Model(1);
Interval.pEquilibrium = test(i).pEquilibrium{1};
Interval.DEquilibrium = test(i).DEquilibrium{1};
Interval.pEfficient = test(i).pEfficient{1};
Interval.DEfficient = test(i).DEfficient{1};
Mandate.Model = test(i).Model(2);
Mandate.pEquilibrium = test(i).pEquilibrium{2};
Mandate.DEquilibrium = test(i).DEquilibrium{2};
Mandate.pEfficient = test(i).pEfficient{2};
Mandate.DEfficient = test(i).DEfficient{2};

nContracts = Interval.Model.nContracts;

xGrid = zeros(1, nContracts);
for j = 1 : nContracts
    xGrid(j) = Interval.Model.meanCoverage(Interval.Model.contracts{j});
end

% Calculate mean loss parameter of agents purchasing each contract
moralHazardTypeVector = zeros(Population(i).size, 1);
meanLossTypeVector = zeros(Population(i).size, 1);
riskAversionVector = zeros(Population(i).size, 1);

for j = 1 : Population(i).size
    moralHazardTypeVector(j) = Population(i).typeList{j}.H;
    meanLossTypeVector(j) = Population(i).typeList{j}.M;
    riskAversionVector(j) = Population(i).typeList{j}.A;
end

[~, ~, ~, choiceVectorEquilibrium] = Population(i).demand(Interval.pEquilibrium);
[~, ~, ~, choiceVectorEfficient]   = Population(i).demand(Interval.pEfficient);
meanLossEquilibrium     = zeros(1, Interval.Model.nContracts);
meanLossEfficient       = zeros(1, Interval.Model.nContracts);

for j = 1 : Interval.Model.nContracts
    buyersEquilibrium = (choiceVectorEquilibrium == j);
    meanLossVectorBuyersEquilibrium = meanLossTypeVector(buyersEquilibrium);
    meanLossEquilibrium(j)  = mean(meanLossVectorBuyersEquilibrium);
    buyersEfficient = choiceVectorEfficient == j;
    meanLossVectorBuyersEfficient = meanLossTypeVector(buyersEfficient);
    meanLossEfficient(j)  = mean(meanLossVectorBuyersEfficient);
end;

slopeVectorEquilibrium = zeros(populationSize, 1);
for j = 1 : nContracts
    I = (choiceVectorEquilibrium == j);
    slopeVectorEquilibrium(I) = Interval.Model.meanCoverage(Interval.Model.contracts{j});
end

save(['figureData' num2str(i)],'-regexp','^(?!(Population|test)$).+')

end



