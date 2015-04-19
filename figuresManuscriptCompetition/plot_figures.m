%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Script to plot figures for the paper. %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Start
clear;
close all;
addpath('../classes');
addpath('../plotFunctions');
addpath('../plotFunctions/export_fig');
addpath('./plotFunctions');
mkdir('figures');
rng(1);

%% Individual models
for modelNameString = { ...
        'interval',                  ...
        'mandate',                   ...
        'interval_high_mh_variance', ...
        'mandate_high_mh_variance'};
    
    CalculationData = load(modelNameString{1});
    plotEquilibriumAndOptimum( ...
        CalculationData.Model, CalculationData.pEquilibrium, CalculationData.DEquilibrium, ...
        CalculationData.pEfficient, CalculationData.DEfficient, modelNameString{1});
end;

%% Interval vs mandate
% Start
    close all;
    clear;
    Interval = load('interval');
    Mandate  = load('mandate');
    
% Calculate necessary series
    nContracts = Interval.Model.nContracts;
    xGrid = zeros(1, Interval.Model.nContracts);
    for j = 1 : nContracts
        xGrid(j) = Interval.Model.contracts{j}.slope;
    end;

% Equilibrium vs mandate prices calculations
    Population = population(Interval.Model, 10^5);
    meanLossTypeVector = zeros(Population.size, 1);
    for i = 1 : Population.size
        meanLossTypeVector(i) = Population.typeList{i}.M;
    end;
    % Prices for the mandate, setting prices for contracts below minimum to
    % infinity.
    pMandate = [zeros(1, nContracts - Mandate.Model.nContracts) + 10^6, Mandate.pEquilibrium];
    DMandate = Population.demand(pMandate);

    % Calculate choice vectors, for the model with all contracts available
    [~, ~, ~, Interval.choiceVector]   = Population.demand(Interval.pEquilibrium);
    [~, ~, ~, Mandate.choiceVector]    = Population.demand(pMandate);
    meanLossEquilibrium     = zeros(1, nContracts);
    meanLossMandate         = zeros(1, nContracts);
    for j = 1 : nContracts
        buyersEquilibrium = (Interval.choiceVector == j);
        meanLossVectorBuyersEquilibrium = meanLossTypeVector(buyersEquilibrium);
        meanLossEquilibrium(j)  = mean(meanLossVectorBuyersEquilibrium);

        buyersMandate = (Mandate.choiceVector == j);
        meanLossVectorBuyersMandate = meanLossTypeVector(buyersMandate);
        meanLossMandate(j)  = mean(meanLossVectorBuyersMandate);
    end;

    % For plotting set prices to Inf below minimum
    pMandate = [zeros(1, nContracts - Mandate.Model.nContracts) + Inf, Mandate.pEquilibrium];
    

% Eqm vs mandate prices graphs
colorRed  = [230,77,79]   / 255;
colorBlue = [208,226,241]   / 255;
colorBlueLine = colorBlue     / 2;
colormap([colorRed; colorBlue]);


    figMandateEqPrices = figure;
        set(figMandateEqPrices, 'Position', [0 0 1.618.*500 500]);
        set(figMandateEqPrices, 'name', 'Price and Mean Loss with and without a mandate', 'numbertitle', 'on');
        % Equilibrium
        % Plot prices
            plot(xGrid, Interval.pEquilibrium, 'k', 'linewidth', 3);
            hold on;
        % Plot average loss parameter
            scatter(xGrid, meanLossEquilibrium, max(1, Interval.DEquilibrium.*3600), ...
                'MarkerFaceColor', colorRed, ...
                'MarkerEdgeColor', 'k', ...
                'LineWidth', 1.5); 
        % Efficient
        % Plot prices
            plot(xGrid, pMandate, 'color', colorBlueLine, 'linewidth', 3);
            hold on;
        % Plot average loss parameter
            scatter(xGrid, meanLossMandate, max(1, DMandate.*3600), ...
                'MarkerFaceColor', colorBlue, ...
                'MarkerEdgeColor', 'k', ...
                'LineWidth', 1.5);        
        % Labels
            legend('No Mandate - Prices', 'No Mandate - Losses', ...
                   'Mandate - Prices', 'Mandate - Losses', ...
                   'Location', 'NorthWest');
            legend boxoff;
            xlabel('Contract');
            ylabel('($)');
       % Other options
            box off;
            set(gca,'FontSize',27);
            set(findall(gcf,'type','text'),'FontSize',27);
       % Financial tick label
            axis([0, 1, 0, 8000])
            set(gca, 'YTickLabel', num2bank(get(gca, 'YTick')));

    fileName = ['./figures/mandate_vs_equilibrium_prices.pdf'];        
    export_fig(fileName, '-transparent');
    fileName = [fileName(1: length(fileName)-4), '.eps'];
    print(fileName, '-depsc2');

% Histograms
    figMandateQuantities = figure;
    set(figMandateQuantities, 'Position', [0 0 1.618.*500 500]);
    set(figMandateQuantities, 'name', 'Quantities with and without a Mandate', 'numbertitle', 'on');
    % Plot quantities
        bar(xGrid, ...
            nContracts .* ...
            [Interval.DEquilibrium; ...
            zeros(1, nContracts - Mandate.Model.nContracts), Mandate.DEquilibrium]', 3);
        colormap([colorRed; colorBlue]);
        hold on;
        axis([0, 1, 0, .1 .* nContracts]);
    % Labels
        legend('No Mandate', 'Mandate', 'Location', 'NorthWest');
        legend boxoff
        xlabel('Coverage');
        ylabel('(Density)');
   % Other options
        box off;
        set(gca,'FontSize',27);
        set(findall(gcf,'type','text'),'FontSize',27);
        
   % Save
       fileName = './figures/mandate_vs_equilibrium_quantities.pdf';        
       export_fig(fileName, '-transparent');
       fileName = [fileName(1: length(fileName)-4), '.eps'];
       print(fileName, '-depsc2');