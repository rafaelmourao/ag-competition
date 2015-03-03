%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_caralognormal_model_nl_utility %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-03-03

This script tests the methods uFunction and cFunction in the caral class.
%}


clear;
close all;
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

% Set up model
deductibleVector  = [1500, 750, 500, 250, 0];
coinsuranceVector = [.1, .1, .1, .1, .1];
oopMaxVector      = [4500, 3750, 3500, 2750, 2500];

Model = healthcaralognormalmodel_nl( ...
    deductibleVector, coinsuranceVector, oopMaxVector, ...
    typeDistributionMean, typeDistributionLogCovariance);

%% Consumer with no moral hazard and small distribution of losses
close all;
x.deductible  = 100;
x.coinsurance = 0.5;
x.oopMax      = 1000;
x.name        = 'test';

type.A = 1e-5;
type.H = 0;
type.M = 3000;
type.S = 10;

display(x);
display(type);

display('Test exPostUtility');
    nGrid = 30;
    lossGrid = linspace(-500, 3000, nGrid);
    uExPost = zeros(1, nGrid);
    eExPost = zeros(1, nGrid);
    for i = 1:nGrid
        [uExPost(i), eExPost(i), limits] = Model.exPostUtility(x, type, lossGrid(i));
    end;
    display(uExPost);
    display(eExPost);
    hold on;
    plot(lossGrid, [eExPost; lossGrid]);
    
    figure();
    plot(lossGrid, [-uExPost]);    

display('Test uFunction and cFunction');
    u    = Model.uFunction(Model.contracts{1}, type);
    c    = Model.cFunction(Model.contracts{1}, type);
    display(u);
    display(c);
    
%% Consumer with no moral hazard and small distribution of losses
close all;
type.H = 1000;

display(x);
display(type);

display('Test exPostUtility');
    nGrid = 20;
    lossGrid = linspace(-500, 3000, nGrid);
    uExPost = zeros(1, nGrid);
    eExPost = zeros(1, nGrid);
    for i = 1:nGrid
        [uExPost(i), eExPost(i), limits] = Model.exPostUtility(x, type, lossGrid(i));
    end;
    display(uExPost);
    display(eExPost);
    hold on;
    plot(lossGrid, [eExPost; lossGrid]);
    
    figure();
    plot(lossGrid, [-uExPost]);    

display('Test uFunction and cFunction');
    u    = Model.uFunction(Model.contracts{1}, type);
    c    = Model.cFunction(Model.contracts{1}, type);
    display(u);
    display(c);    