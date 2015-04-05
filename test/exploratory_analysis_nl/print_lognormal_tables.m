%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print_lognormal_tables.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-03-20

Simple script to print test tables from variables. Create matlab tables and
print to text files

%}

addpath('../../classes')
load('tests.mat')
format shortg

file = @(x) ['test_lognormal_' num2str(x) '.txt'];

for i = 1:length(test)
   
ncontracts = length(test(i).contracts);
nPopulations = test(i).nPopulations;
contractinfo = {'Deductible','Coinsurance','OOP Max','Mean Coverage'}';
for j = 1:ncontracts
    contractinfo{1,j+1} = test(i).Model(1).contracts{j}.deductible;
    contractinfo{2,j+1} = test(i).Model(1).contracts{j}.coinsurance; 
    contractinfo{3,j+1} = test(i).Model(1).contracts{j}.oopMax;
    contractinfo{4,j+1} = round( test(i).Model(1).meanCoverage(test(i).Model(1).contracts{j}), 2);
end

fid = fopen(file(i),'w+');
string = evalc('disp(contractinfo)');
fwrite(fid,string);

fprintf(fid,'\nEquilibrium Calculation Parameters\n\n');
fwrite(fid,evalc('disp(test(i).CalculationParametersEquilibrium)'));
fprintf(fid,'\nOptimum Calculation Parameters\n\n');
fwrite(fid,evalc('disp(test(i).CalculationParametersOptimum)'));

x = cell(8,ncontracts);
Tables = cell(1,ncontracts);


for j = 1:nPopulations
    x(1,:)=num2cell(round(test(i).DEquilibrium{j},2));
    x(2,:)=num2cell(test(i).pEquilibrium{j});
    x(3,1)=num2cell(test(i).WEquilibrium{j});
    x(4,1)=num2cell(test(i).ComputationOutputEquilibrium{j}.error);
    x(5,:)=num2cell(round(test(i).DEfficient{j},2));
    x(6,:)=num2cell(test(i).pEfficient{j});
    x(7,1)=num2cell(test(i).WEfficient{j});
    x(8,1)=num2cell(test(i).ComputationOutputEfficient{j}.error);
    x(9,test(i).mandate:end)=num2cell(round(test(i).DEquilibrium{j+nPopulations},2));
    x(10,test(i).mandate:end)=num2cell(test(i).pEquilibrium{j+nPopulations});
    x(11,1)=num2cell(test(i).WEquilibrium{j+nPopulations});
    x(12,1)=num2cell(test(i).ComputationOutputEfficient{j+nPopulations}.error);

    rowtitles = {'Equilibrium Demand','Equilibrium Prices',...
          'Equilibrium Welfare', 'Equilibrium Error', 'Efficient Demand','Efficient Prices',...
          'Efficient Welfare','Efficient Error', 'Mandate Demand','Mandate Prices','Mandate Welfare','Mandate Error'}';
      
    string = evalc('disp([rowtitles,x])');
    
%     fprintf(fid,['\n\n' Testnames{j} '\n\n\n']);
    fprintf(fid,'\n\nMean Parameters\n\n');
    fwrite(fid,evalc('disp(num2cell(test(i).Model(j).typeDistributionMean))'));
    fprintf(fid,'\nLog Covariances\n\n');
    fwrite(fid,evalc('disp(num2cell(test(i).Model(j).typeDistributionLogCovariance))'));
    fwrite(fid,string);
end

fclose(fid);

end

