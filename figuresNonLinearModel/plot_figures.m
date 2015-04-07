clear;
close all;
mkdir('figures');
addpath('../classes');
addpath('../plotFunctions');
addpath('../plotFunctions/export_fig');

for i = 1:3

 close all
 load(['figureData' num2str(i) '.mat'])

% Parameters
colorRed  = [230,77,79]   / 255;
colorBlue = [208,226,241]   / 255;
colorBlueLine = colorBlue     / 2;

% Equilibrium and optimum
figEqPrices = figure;
    set(figEqPrices, 'Position', [0 0 1.618.*500 500]);
    set(figEqPrices, 'name', 'Equilibrium Price and Mean Loss', 'numbertitle', 'on');
    % Plot prices
        plot(xGrid, Interval.pEquilibrium, 'k', 'linewidth', 3);
        hold on;
    % Plot average loss parameter
        scatter(xGrid, meanLossEquilibrium, max(1,Interval.DEquilibrium.*3600)*0.85, ...
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

fileName = ['./figures/', num2str(Interval.Model.typeDistributionLogCovariance(1)), '_', 'equilibrium_prices.pdf'];        
export_fig(fileName, '-transparent');
fileName = [fileName(1: length(fileName)-4), '.eps'];
print(fileName, '-depsc2');
        

figEqDemandProfile = figure;
    set(figEqDemandProfile, 'Position', [0 0 1.618.*500 500]);
    set(figEqDemandProfile, 'name', 'Equilibrium Demand Profile', 'numbertitle', 'on');
    % Plot scatter
        scatter(moralHazardTypeVector(1:nPopulationDemandProfile), ...
            riskAversionVector(1:nPopulationDemandProfile), 20, ...
            slopeVectorEquilibrium(1:nPopulationDemandProfile));
        colorbar;
            caxis([0.01, 1]);
    % Labels
        xlabel('Moral Hazard,  H_{\theta}');
        ylabel('Risk Aversion, A_{\theta}');
    % Properties
        set(gca, 'xscale', 'log', ...
                 'yscale', 'log');
        grid off;
   % Other options
        box off;
        set(gca,'FontSize',27);
        set(findall(gcf,'type','text'),'FontSize',27);        
    % Financial tick label
        ylim([10^(-5),10^(-1)]) 

fileName = ['./figures/', num2str(Interval.Model.typeDistributionLogCovariance(1)), '_', 'scatter_a_h.pdf'];        
export_fig(fileName, '-transparent');
fileName = [fileName(1: length(fileName)-4), '.eps'];
print(fileName, '-depsc2');

figMandateQuantities = figure;
set(figMandateQuantities, 'Position', [0 0 1.618.*500 500]);
set(figMandateQuantities, 'name', 'Quantities with and without a Mandate', 'numbertitle', 'on');
% Plot quantities
bar(xGrid, ...
           [Interval.DEquilibrium; ...
            zeros(1, nContracts - Mandate.Model.nContracts), Mandate.DEquilibrium]');
        colormap([colorRed; colorBlue]);
        hold on;
        axis([min(xGrid), 1, 0, 1]);
    % Labels
        legend('No Mandate', 'Mandate', 'Location', 'NorthWest');
        legend boxoff
        xlabel('Coverage');
        ylabel('Frequency');
   % Other options
        box off;
        set(gca,'FontSize',27);
        set(findall(gcf,'type','text'),'FontSize',27);
   % Save
       fileName = ['./figures/' num2str(Interval.Model.typeDistributionLogCovariance(1)) '_mandate_vs_equilibrium_quantities.pdf'];        
       export_fig(fileName, '-transparent');
       fileName = [fileName(1: length(fileName)-4), '.eps'];
print(fileName, '-depsc2');

end