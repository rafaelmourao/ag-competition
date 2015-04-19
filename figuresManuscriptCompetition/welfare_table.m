clear;
addpath('../classes');
rng(1);

% Draws a table with welfare comparisons. Saves a .csv for high moral
% hazard and a .csv for low moral hazard.

% Table structure:

%          Eqm p      Optimal p
% X      W   E[x]     W   E[x]
%        % Low MH %
% [0,1]
% [60,1]
% [0 90]
% [0 60 90]
% [60 90]
%        % High MH %
% [0,1]
% [60,1]
% [0 90]
% [0 60 90]
% [60 90]

tableLowMH  = zeros(5, 4);
tableHighMH = zeros(5, 4);

% [0, 1] Low
Interval = load('interval');
slopeVector = zeros(1, Interval.Model.nContracts);
for ii = 1 : Interval.Model.nContracts
    slopeVector(ii) = Interval.Model.contracts{ii}.slope;
end;
tableLowMH(1, :) = [ ...
    Interval.WEquilibrium, ...
    slopeVector * Interval.DEquilibrium', ...
    Interval.WEfficient, ...
    slopeVector * Interval.DEfficient'];

% [0, 1] High
IntervalHighMH = load('interval_high_mh_variance');
slopeVector = zeros(1, IntervalHighMH.Model.nContracts);
for ii = 1 : IntervalHighMH.Model.nContracts
    slopeVector(ii) = IntervalHighMH.Model.contracts{ii}.slope;
end;
tableHighMH(1, :) = [ ...
    IntervalHighMH.WEquilibrium, ...
    slopeVector * IntervalHighMH.DEquilibrium', ...
    IntervalHighMH.WEfficient, ...
    slopeVector * IntervalHighMH.DEfficient'];

% [0.60 1] Low
Mandate = load('mandate');
slopeVector = zeros(1, Mandate.Model.nContracts);
for ii = 1 : Mandate.Model.nContracts
    slopeVector(ii) = Mandate.Model.contracts{ii}.slope;
end;
tableLowMH(2, :) = [ ...
    Mandate.WEquilibrium, ...
    slopeVector * Mandate.DEquilibrium', ...
    Mandate.WEfficient, ...
    slopeVector * Mandate.DEfficient'];

% [0.60 1] High
MandateHighMH = load('mandate_high_mh_variance');
slopeVector = zeros(1, MandateHighMH.Model.nContracts);
for ii = 1 : MandateHighMH.Model.nContracts
    slopeVector(ii) = MandateHighMH.Model.contracts{ii}.slope;
end;
tableHighMH(2, :) = [ ...
    MandateHighMH.WEquilibrium, ...
    slopeVector * MandateHighMH.DEquilibrium', ...
    MandateHighMH.WEfficient, ...
    slopeVector * MandateHighMH.DEfficient'];

% 60 90
for tableRow = [3 4 5]
    if tableRow == 3
        slopeVector = [0.00 0.90];
    elseif tableRow == 4
        slopeVector = [0.60 0.90];
    elseif tableRow == 5
        slopeVector = [0 0.60 0.90];
    end;
    

    % Low MH variance
    Model = healthcaralognormalmodel(slopeVector, ...
        Interval.Model.typeDistributionMean, Interval.Model.typeDistributionLogCovariance);
    rng(1);
    Population = population(Model, Interval.populationSize);

    [pEquilibrium, DEquilibrium, ACEquilibrium, ComputationOutputEquilibrium] = ...
        Population.findequilibrium(Interval.CalculationParametersEquilibrium);
    WEquilibrium = Population.welfare(pEquilibrium, Interval.costOfPublicFunds);

    [pEfficient, WEfficient, ComputationOutputEfficient] = ...
        findefficient(Population, Interval.costOfPublicFunds, Interval.CalculationParametersOptimum);
    DEfficient = Population.demand(pEfficient);

    tableLowMH(tableRow, :) = [ ...
        WEquilibrium, ...
        slopeVector * DEquilibrium', ...
        WEfficient, ...
        slopeVector * DEfficient'];
    
    % High MH variance
    Model = healthcaralognormalmodel(slopeVector, ...
        IntervalHighMH.Model.typeDistributionMean, IntervalHighMH.Model.typeDistributionLogCovariance);
    Population = population(Model, Interval.populationSize);

    [pEquilibrium, DEquilibrium, ACEquilibrium, ComputationOutputEquilibrium] = ...
        Population.findequilibrium(Interval.CalculationParametersEquilibrium);
    WEquilibrium = Population.welfare(pEquilibrium, Interval.costOfPublicFunds);

    [pEfficient, WEfficient, ComputationOutputEfficient] = ...
        findefficient(Population, Interval.costOfPublicFunds, Interval.CalculationParametersOptimum);
    DEfficient = Population.demand(pEfficient);

    tableHighMH(tableRow, :) = [ ...
        WEquilibrium, ...
        slopeVector * DEquilibrium', ...
        WEfficient, ...
        slopeVector * DEfficient'];
end;

% Normalize welfare relative to Laissez-faire
tableLowMH(:, [1, 3])  = tableLowMH(:, [1, 3])  - tableLowMH(1, 1);
tableHighMH(:, [1, 3]) = tableHighMH(:, [1, 3]) - tableHighMH(1, 1);

% Save latex table.
latexHeader = '\begin{tabular}{lccccccccccc}     & \multicolumn{5}{c}{$\sigma^2_H=0.28$}                              &  & \multicolumn{5}{c}{$\sigma^2_H=0.98$}                              \\ \cline{2-6} \cline{8-12}      & \multicolumn{2}{c}{Equilibrium} &  & \multicolumn{2}{c}{Efficient} &  & \multicolumn{2}{c}{Equilibrium} &  & \multicolumn{2}{c}{Efficient} \\ $X$ & Welfare         & $E[x]$        &  & Welfare        & $E[x]$       &  & Welfare         & $E[x]$        &  & Welfare        & $E[x]$      ';
latexEnd    = '\end{tabular}';

matrixColumn = -99 * ones(5, 1);
matrixTable  = [matrixColumn, ...
    tableLowMH(:, 1:2), matrixColumn, tableLowMH(:, 3:4), matrixColumn, ...
    tableHighMH(:, 1:2), matrixColumn, tableHighMH(:, 3:4)]; 
cellTable = num2cell(matrixTable);

latexMiddle = '';
for iLine = 1:5
    for jColumn = 1:12
        if cellTable{iLine, jColumn} == -99
            cellTable{iLine, jColumn} = '';
        end;
        
        if (jColumn ~= 2) && (jColumn ~= 5) && (jColumn ~= 8) && (jColumn ~= 11)
            cellTable{iLine, jColumn} = num2str(cellTable{iLine, jColumn}, '%.2f');
        else
            cellTable{iLine, jColumn} = num2str(cellTable{iLine, jColumn}, '%.0f');
        end;
        
        if jColumn == 1
            switch iLine
                case 1
                    cellTable{iLine, jColumn} = '$[0,1]$';
                case 2
                    cellTable{iLine, jColumn} = '$[0.60,1]$';
                case 3
                    cellTable{iLine, jColumn} = '0, 0.90';
                case 4
                    cellTable{iLine, jColumn} = '0.60, 0.90';
                case 5
                    cellTable{iLine, jColumn} = '0, 0.60, 0.90';
            end;
        end;
        
        if jColumn == 1
            latexMiddle = [latexMiddle, ' \\ ', cellTable{iLine, jColumn}];
        else
            latexMiddle = [latexMiddle, ' & ', cellTable{iLine, jColumn}];
        end;
    end;
end;

latexCode = [latexHeader, latexMiddle, latexEnd];

% Save
fName = 'welfare_table.tex';
fId   = fopen(fName, 'w');
fprintf(fId, '%s\r\n', latexCode);
fclose(fId);

clear Population;
save('welfare_table');

% Save welfare gain numbers
welfareGainMandate = tableLowMH(2, 1);
welfareGainOptimal = tableLowMH(1, 3);

fileID = fopen('welfare_table_gain_mandate.tex', 'w');
fprintf(fileID, '%0.0f', welfareGainMandate);
fclose(fileID);

fileID = fopen('welfare_table_gain_optimal_regulation.tex', 'w');
fprintf(fileID, '%0.0f', welfareGainOptimal);
fclose(fileID);