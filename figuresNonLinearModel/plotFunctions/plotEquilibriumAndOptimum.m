function [  ] = plotEquilibriumAndOptimum( ...
    Model, pEquilibrium, DEquilibrium, pEfficient, DEfficient, modelNameString, Population)
%plotEquilibriumAndOptimum Plots graphs that we sometimes use in every model.
%   Inputs are the model, equilibrium prices and demand, efficient prices
%   and demand, and a string with the name of the model. This function
%   generates no outputs. Instead it plots equilibrium prices with mean
%   losses, same for efficient, a graph overlaying efficient and
%   equilibrium, and a demand profile. You have to change in the function
%   what population size to use for plots. I typically set it to 10000.

% Start
rng(1);

% Parameters
if nargin < 7
populationSize           = 10^5;
nPopulationDemandProfile = min(7000, populationSize);
Population = population(Model, populationSize);
else
populationSize           = Population.size;
nPopulationDemandProfile = min(7000, populationSize);
end

% Color Scheme
colorRed  = [230,77,79]   / 255;
colorBlue = [208,226,241]   / 255;
colorBlueLine = colorBlue     / 2;

% Prepare
xGrid = zeros(1, Model.nContracts);
for j = 1 : Model.nContracts
    xGrid(j) = Model.meanCoverage(Model.contracts{j});
end;

% Calculate mean loss parameter of agents purchasing each contract
meanLossTypeVector = zeros(Population.size, 1);
riskAversionVector = zeros(Population.size, 1);
for i = 1 : Population.size
    meanLossTypeVector(i) = Population.typeList{i}.M;
    riskAversionVector(i) = Population.typeList{i}.A;
end;

[~, ~, ~, choiceVectorEquilibrium] = Population.demand(pEquilibrium);
[~, ~, ~, choiceVectorEfficient]   = Population.demand(pEfficient);
meanLossEquilibrium     = zeros(1, Model.nContracts);
meanLossEfficient       = zeros(1, Model.nContracts);

for j = 1 : Model.nContracts
    buyersEquilibrium = (choiceVectorEquilibrium == j);
    meanLossVectorBuyersEquilibrium = meanLossTypeVector(buyersEquilibrium);
    meanLossEquilibrium(j)  = mean(meanLossVectorBuyersEquilibrium);
    
    buyersEfficient = choiceVectorEfficient == j;
    meanLossVectorBuyersEfficient = meanLossTypeVector(buyersEfficient);
    meanLossEfficient(j)  = mean(meanLossVectorBuyersEfficient);
end;

slopeVectorEquilibrium = zeros(populationSize, 1);
for j = 1 : Model.nContracts
    I = (choiceVectorEquilibrium == j);
    slopeVectorEquilibrium(I) = Model.meanCoverage(Model.contracts{j});
end;

% Equilibrium and optimum
figEqPrices = figure;
    set(figEqPrices, 'Position', [0 0 1.618.*500 500]);
    set(figEqPrices, 'name', 'Equilibrium Price and Mean Loss', 'numbertitle', 'on');
    % Plot prices
        plot(xGrid, pEquilibrium, 'k', 'linewidth', 3);
        hold on;
    % Plot average loss parameter
        scatter(xGrid, meanLossEquilibrium, max(1,DEquilibrium.*3600), ...
            'MarkerFaceColor', colorRed, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.5);
    % Labels
        legend('Equilibrium Prices', 'Average Loss Parameter', 'Location', 'NorthWest');
        legend boxoff
        xlabel('Contract');
        ylabel('($)');
   % Other options
        box off;
        set(gca,'FontSize',27);
        set(findall(gcf,'type','text'),'FontSize',27);
   % Financial tick label
        axis([min(xGrid),1,0,8000])
        set(gca, 'YTickLabel', num2bank(get(gca, 'YTick')));

fileName = ['./figures/', modelNameString, '_', 'equilibrium_prices.pdf'];        
export_fig(fileName, '-transparent');
fileName = [fileName(1: length(fileName)-4), '.eps'];
print(fileName, '-depsc2');
        
        
figEffPrices = figure;
    set(figEffPrices, 'Position', [0 0 1.618.*500 500]);
    set(figEffPrices, 'name', 'Optimum Price and Mean Loss', 'numbertitle', 'on');
    % Plot prices
        plot(xGrid, pEfficient, 'k', 'linewidth', 3);
        hold on;
    % Plot average loss parameter
        scatter(xGrid, meanLossEfficient, max(1,DEfficient.*3600), ...
            'MarkerFaceColor', colorRed, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.5);
    % Labels
        legend('Optimum Prices', 'Average Loss Parameter', 'Location', 'NorthWest');
        legend boxoff;
        xlabel('Contract');
        ylabel('($)');
   % Other options
        box off;
        set(gca,'FontSize',27);
        set(findall(gcf,'type','text'),'FontSize',27);
   % Financial tick label
        axis([min(xGrid),1,0,8000])
        set(gca, 'YTickLabel', num2bank(get(gca, 'YTick')));  

fileName = ['./figures/', modelNameString, '_', 'efficient_prices.pdf'];        
export_fig(fileName, '-transparent');
fileName = [fileName(1: length(fileName)-4), '.eps'];
print(fileName, '-depsc2');

figEffnEqPrices = figure;
    set(figEffnEqPrices, 'Position', [0 0 1.618.*500 500]);
    set(figEffnEqPrices, 'name', 'Optimum and Equilibrium Price and Mean Loss', 'numbertitle', 'on');
    % Equilibrium
    % Plot prices
        plot(xGrid, pEquilibrium, 'k', 'linewidth', 3);
        hold on;
    % Plot average loss parameter
        scatter(xGrid, meanLossEquilibrium, max(1,DEquilibrium.*3600), ...
            'MarkerFaceColor', colorRed, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.5); 
    % Efficient
    % Plot prices
        plot(xGrid, pEfficient, 'color', colorBlueLine, 'linewidth', 3);
        hold on;
    % Plot average loss parameter
        scatter(xGrid, meanLossEfficient, max(1,DEfficient.*3600), ...
            'MarkerFaceColor', colorBlue, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.5);        
    % Labels
        legend('Equilibrium Prices', 'Equilibrium Losses', ...
               'Optimum Prices', 'Optimum Losses', ...
               'Location', 'NorthWest');
        legend boxoff;
        xlabel('Contract');
        ylabel('($)');
   % Other options
        box off;
        set(gca,'FontSize',27);
        set(findall(gcf,'type','text'),'FontSize',27);
   % Financial tick label
        axis([min(xGrid),1,0,8000])
        set(gca, 'YTickLabel', num2bank(get(gca, 'YTick')));
        
fileName = ['./figures/', modelNameString, '_', 'efficient_and_equilibrium_prices.pdf'];        
export_fig(fileName, '-transparent');
fileName = [fileName(1: length(fileName)-4), '.eps'];
print(fileName, '-depsc2');

end