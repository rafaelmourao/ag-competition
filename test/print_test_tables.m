%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print_test_tables.m      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-03-10

Simple script to print test tables from variables. Create matlab tables and
print to an Excel file

%}

addpath('../classes')
file = @(x) ['test' num2str(x)];
x = cell(6,6);
Tables = cell(1,6);

Testnames = {'High Aversion, Low MH Variance', 'Low Aversion, Low MH Variance' ...
    'High Aversion, High MH Variance', 'Low Aversion, High MH Variance'};
for i = 1:4
    load(file(i));
    for j = 1:6;
    a(j) = Model.meanCoverage(Model.contracts{j});
    end
    x(1,:)=num2cell(DEquilibrium);
    x(2,:)=num2cell(pEquilibrium);
    x(3,:)=num2cell(DEfficient);
    x(4,:)=num2cell(pEfficient);
    x(5,1)=num2cell(WEfficient);
    load(file(i+4))
    x(6,4:end)=num2cell(DEquilibrium);
    x(7,4:end)=num2cell(pEquilibrium);
    x(8,:)=num2cell(a);
    
    
    Tables{i} = cell2table(x,...
          'RowNames',{'Equilibrium Demand','Equilibrium Prices',...
         'Efficient Demand','Efficient Prices','Efficient Welfare',...
         'Mandate Demand','Mandate Prices','Mean Coverage'},...
         'VariableNames',strcat('Contract ',{'0','1','2','3','4','5'}));
     
     writetable(Tables{i},['Tables.xls'],'Sheet',Testnames{i},...
         'WriteRowNames',true)
end


