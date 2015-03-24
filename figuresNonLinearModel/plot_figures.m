%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Script to plot figures for the paper. %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Start
clear;
close all;
mkdir('figures');
addpath('../classes');
addpath('../plotFunctions');
addpath('../plotFunctions/export_fig');
addpath('./plotFunctions');
rng(1);

%% Interval vs mandate


% Start
    load('Populations.mat')
    load('tests.mat')
    Interval.Model = test.Model(1);
    Interval.pEquilibrium = test.pEquilibrium{1};
    Interval.DEquilibrium = test.DEquilibrium{1};
    Interval.pEfficient = test.pEfficient{1};
    Interval.DEfficient = test.DEfficient{1};
    Mandate.Model = test.Model(2);
    Mandate.pEquilibrium = test.pEquilibrium{2};
    Mandate.DEquilibrium = test.DEquilibrium{2};
    Mandate.pEfficient = test.pEfficient{2};
    Mandate.DEfficient = test.DEfficient{2};
    
    plotEquilibriumAndOptimum( ...
    Interval.Model, Interval.pEquilibrium, Interval.DEquilibrium, ...
    Interval.pEfficient, Interval.DEfficient, 'Interval', Population);

    % Calculate necessary series

    nContracts = Interval.Model.nContracts;
    xGrid = zeros(1, Interval.Model.nContracts);
    for j = 1 : nContracts
        xGrid(j) = Interval.Model.meanCoverage(Interval.Model.contracts{j});
    end

    meanLossTypeVector = zeros(Population.size, 1);
    for i = 1 : Population.size
        
        meanLossTypeVector(i) = Population.typeList{i}.M;
    end
    
    % Prices for the mandate, setting prices for contracts below minimum to
    % infinity.
    pMandate = [zeros(1, nContracts - Mandate.Model.nContracts) + 10^6, test.pEquilibrium{2}];
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
    pMandate = [zeros(1, nContracts - test.Model(2).nContracts) + Inf, Mandate.pEquilibrium];
    

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
            axis([min(xGrid), 1, 0, 8000])
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
        axis([min(xGrid), 1, 0, .1 .* nContracts]);
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