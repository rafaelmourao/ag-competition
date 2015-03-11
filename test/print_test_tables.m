%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print_test_tables.m      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-03-10

Simple script to print test tables from variables. Create matlab tables and
print to an Excel file

%}


file = @(x) ['test' num2str(x)];
x = cell(6,6);
Tables = cell(1,6);

Testnames = {'High Aversion, Low MH Variance', 'Low Aversion, Low MH Variance' ...
    'High Aversion, High MH Variance', 'Low Aversion, High MH Variance'};
for i = 1:4
    load(file(i));
    x(1,:)=num2cell(DEquilibrium);
    x(2,:)=num2cell(pEquilibrium);
    x(3,:)=num2cell(DEfficient);
    x(4,:)=num2cell(pEfficient);
    load(file(i+4))
    x(5,4:end)=num2cell(DEquilibrium);
    x(6,4:end)=num2cell(pEquilibrium);
    
    Tables{i} = cell2table(x,...
          'RowNames',{'Equilibrium Demand','Equilibrium Prices',...
         'Efficient Demand','Efficient Prices',...
         'Mandate Demand','Mandate Prices'},...
         'VariableNames',strcat('Contract ',{'1','2','3','4','5','6'}));
     
     writetable(Tables{i},['Tables.xls'],'Sheet',Testnames{i},...
         'WriteRowNames',true)
end


