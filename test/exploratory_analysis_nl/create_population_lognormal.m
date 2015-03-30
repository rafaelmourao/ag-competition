%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create_population_lognormal.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-03-18
%}

clear;
addpath('../../classes')
rng(1);

nworkers = 128;

% Input model parameters
meanS = sqrt(25000^2 - 5100^2);
typeDistributionMean{1} = ...
    [1.9*1e-3, 1330, 4340, meanS]
typeDistributionMean{2} = ...
    [1.9*1e-4, 1330, 4340, meanS]
typeDistributionMean{3} = ...
    [1e-4, 1330, 4340, meanS]
typeDistributionMean{4} = ...
    [1.9*1e-3, 1330, 4340, meanS]
typeDistributionMean{5} = ...
    [1.9*1e-4, 1330, 4340, meanS]
typeDistributionMean{6} = ...
    [1e-4, 1330, 4340, meanS]
typeDistributionMean{7} = ...
    [1.9*1e-3, 1330, 4340, meanS]
typeDistributionMean{8} = ...
    [1.9*1e-4, 1330, 4340, meanS]
typeDistributionMean{9} = ...
    [1e-4, 1330, 4340, meanS]
typeDistributionMean{10} = ...
    [1.9*1e-3, 1330, 4340, meanS]
typeDistributionMean{11} = ...
    [1.9*1e-4, 1330, 4340, meanS]
typeDistributionMean{12} = ...
    [1e-4, 1330, 4340, meanS]
typeDistributionMean{13} = ...
    [1.9*1e-3, 1330, 4340, meanS]
typeDistributionMean{14} = ...
    [1.9*1e-4, 1330, 4340, meanS]
typeDistributionMean{15} = ...
    [1e-4, 1330, 4340, meanS]
typeDistributionMean{16} = ...
    [1.9*1e-3, 1330, 4340, meanS]
typeDistributionMean{17} = ...
    [1.9*1e-4, 1330, 4340, meanS]
typeDistributionMean{18} = ...
    [1e-4, 1330, 4340, meanS]
typeDistributionMean{19} = ...
    [1.9*1e-3, 1330, 4340, meanS]
typeDistributionMean{20} = ...
    [1.9*1e-4, 1330, 4340, meanS]
typeDistributionMean{21} = ...
    [1e-4, 1330, 4340, meanS]


typeDistributionLogCovariance{1} = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{2} = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{3} = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{4} = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.98 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{5} = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.98 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{6} = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.98 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{7} = ...
    [ 1 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{8} = ...
    [ 1 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{9} = ...
    [ 1 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{10} = ...
    [ 1 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.98 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{11} = ...
    [ 1 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.98 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{12} = ...
    [ 1 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.98 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{13} = ...
    [ 2 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{14} = ...
    [ 2 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{15} = ...
    [ 2 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{16} = ...
    [ 2 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.98 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{17} = ...
    [ 2 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.98 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{18} = ...
    [ 2 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.98 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{19} = ...
    [ 1 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.01 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{20} = ...
    [ 1 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.01 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

typeDistributionLogCovariance{21} = ...
    [ 1 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.01 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25] % ???

costOfPublicFunds = 0;

% Calculation parameters
populationSize = 1e5
publicInsuranceMaximum = 30000

% Contracts

deductibleVector = [ 30000 23000 17000 12000 9200 6700 5200 3300 2100 1300 600 300 0]
coinsuranceVector = [ 1 0.82  0.65  0.52  0.46  0.39 0.35  0.26  0.2  0.16  0.13  0.12 0.1]
oopMaxVector = [ 30000 23500 17500 13000 10000 7800 6400  5000  4200  3500  2900  2700 2500]

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

