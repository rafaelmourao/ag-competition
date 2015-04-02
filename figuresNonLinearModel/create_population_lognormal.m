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
    [1.9*1e-3, 1330, 4340, meanS];
typeDistributionMean{2} = ...
    [1.9*1e-3, 1330, 4340, meanS];
typeDistributionMean{3} = ...
    [1.9*1e-3, 1330, 4340, meanS];


typeDistributionLogCovariance{1} = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.98 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25]; % ???

typeDistributionLogCovariance{2} = ...
    [ 1 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.98 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25]; % ???

typeDistributionLogCovariance{3} = ...
    [ 2 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.98 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25]; % ???

costOfPublicFunds = 0;

% Calculation parameters
populationSize = 1e6;
publicInsuranceMaximum = 30000;

% Contracts

deductibleVector = [ 30000 23000 17000 12000 9200 6700 5200 3300 2100 1300 600 300 0];
coinsuranceVector = [ 1 0.82  0.65  0.52  0.46  0.39 0.35  0.26  0.2  0.16  0.13  0.12 0.1];
oopMaxVector = [ 30000 23500 17500 13000 10000 7800 6400  5000  4200  3500  2900  2700 2500];

% Calculate equilibrium and optimal prices.

if ~isempty(gcp('nocreate'))
    delete(gcp)
end

poolobj = parpool(nworkers);

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
    fprintf('Time to create sample %d: %f\n',i,time)
    
end

delete(poolobj)
clear poolobj

popsize=whos('Population');
popsize = popsize.bytes;
if popsize < 2e9
    save('Populations.mat')
else
    save('Populations.mat','-v7.3')
end

