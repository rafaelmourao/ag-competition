%% Start
% Load data
clear;
close all;
addpath('../agSim');
rng(1);

Mandate_60 = load('mandate');
Mandate_64 = load('mandate_64.mat');
Mandate_40 = load('mandate_40.mat');
Mandate_44 = load('mandate_44.mat');

% Parameters
populationSize = 10^5;

% Constants
spacing = 0.04;

%% Hack for model structure to work properly
Model_aux_60 = healthcaralognormalmodel(0.6:spacing:1, ...
    Mandate_60.Model.typeDistributionMean, ...
    Mandate_60.Model.typeDistributionLogCovariance);
Model_aux_40 = healthcaralognormalmodel(0.4:spacing:1, ...
    Mandate_40.Model.typeDistributionMean, ...
    Mandate_40.Model.typeDistributionLogCovariance);

%% Calculations
% Create population
Population_60 = population(Model_aux_60, populationSize);
Population_40 = population(Model_aux_40, populationSize);

% Calculate vectors of choices
[~, ~, ~, choice_mandate_60] = demand(Population_60, Mandate_60.pEquilibrium);
[~, ~, ~, choice_mandate_40] = demand(Population_40, Mandate_40.pEquilibrium);

% Calculate change in prices
deltaP_60 = (Mandate_64.pEquilibrium(1) - Mandate_60.pEquilibrium(2))/spacing;
deltaP_40 = (Mandate_44.pEquilibrium(1) - Mandate_40.pEquilibrium(2))/spacing;

% Calculate marginal costs
% 60
I = find(choice_mandate_60 == 1);
mc = zeros(length(I), 1);
for ii = I
    type = Population_60.typeList{ii};
    mc(ii) = type.M + 2 * 0.6 * type.H;
end;

mean_mc_60 = mean(mc);

% 40
I = find(choice_mandate_40 == 1);
mc = zeros(length(I), 1);
for ii = I
    type = Population_60.typeList{ii};
    mc(ii) = type.M + 2 * 0.4 * type.H;
end;

mean_mc_40 = mean(mc);

% Comparative statics
% 60
dp_m = (Mandate_60.pEquilibrium(2) - Mandate_60.pEquilibrium(1)) / spacing;
S_I_m_60 = dp_m - mean_mc_60;

% 40
dp_m = (Mandate_40.pEquilibrium(2) - Mandate_40.pEquilibrium(1)) / spacing;
S_I_m_40 = dp_m - mean_mc_40;


%% Save constants
% C_S_I_m_60
    fileID = fopen('C_S_I_m.tex', 'w');
    fprintf(fileID, '%0.0f', S_I_m_60);
    fclose(fileID);

% C_S_I_m_over_100_60
    fileID = fopen('C_S_I_m_over_100.tex', 'w');
    fprintf(fileID, '%0.0f', S_I_m_60/100);
    fclose(fileID);

% C_change_in_prices_actual_60
    fileID = fopen('C_change_in_prices_actual.tex', 'w');
    fprintf(fileID, '%0.0f', -deltaP_60(1) * spacing);
    fclose(fileID);
    
% C_change_in_prices_predicted_60
    fileID = fopen('C_change_in_prices_predicted.tex', 'w');
    fprintf(fileID, '%0.0f', S_I_m_60 * spacing);
    fclose(fileID);

% C_mass_minimum_coverage_40
    fileID = fopen('C_mass_minimum_coverage_40.tex', 'w');
    fprintf(fileID, '%0.0f', 100* Mandate_40.DEquilibrium(1));
    fclose(fileID);

% C_change_in_prices_predicted_40
    fileID = fopen('C_change_in_prices_predicted_40.tex', 'w');
    fprintf(fileID, '%0.0f', S_I_m_40 * spacing);
    fclose(fileID);

% C_change_in_prices_actual_40
    fileID = fopen('C_change_in_prices_actual_40.tex', 'w');
    fprintf(fileID, '%0.0f', -deltaP_40(1) * spacing);
    fclose(fileID);