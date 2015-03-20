%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print_lognormal_tables.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-03-20

Simple script to print test tables from variables. Create matlab tables and
print to an Excel file

%}

addpath('../classes')
load('tests.mat')

file = @(x) ['test_lognormal_' num2str(x) '.xls'];
Testnames = {'High Aversion, Low MH Variance', 'Low Aversion, Low MH Variance' ...
    'High Aversion, High MH Variance', 'Low Aversion, High MH Variance'};

for i = 1:length(test)
    
Tables = cell(1,6);
ncontracts = test(i).Population(1).nContracts;
x = cell(9,ncontracts);
 
for j = 1:4
    a=[];
    for z = 1:ncontracts;
    a(z) = test(i).Model(j).meanCoverage(test(i).Model(j).contracts{z});
    end
    x(1,:)=num2cell(test(i).DEquilibrium{j});
    x(2,:)=num2cell(test(i).pEquilibrium{j});
    x(3,1)=num2cell(test(i).WEquilibrium{j});
    x(4,:)=num2cell(test(i).DEfficient{j});
    x(5,:)=num2cell(test(i).pEfficient{j});
    x(6,1)=num2cell(test(i).WEfficient{j});
    x(7,2:end)=num2cell(test(i).DEquilibrium{j+4});
    x(8,2:end)=num2cell(test(i).pEquilibrium{j+4});
    x(9,:)=num2cell(a);
    
    
    Tables{j} = cell2table(x,...
          'RowNames',{'Equilibrium Demand','Equilibrium Prices',...
         'Equilibrium Welfare', 'Efficient Demand','Efficient Prices',
         'Efficient Welfare','Mandate Demand','Mandate Prices','Mean Coverage'});
     
     writetable(Tables{j},file(i),'Sheet',Testnames{j},...
         'WriteRowNames',true)
end

end

