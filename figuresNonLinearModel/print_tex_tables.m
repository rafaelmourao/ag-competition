%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print_tex_tables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-04-04

Simple script to print tex tables from variables. Create tex files with
contracts and results tables. 

%}

addpath('../classes')
mkdir('tables')
load('tests.mat')
contract_titles = {'0 (public insurance)','1','2','3','4','5','6 (bronze)','7','8',...
    '9','10','11','12 (Einav et al.)'};

%% Writing contracts table
fid = fopen('tables/contracts.tex','w+');
x = cell(13,5);
x(:,1) = contract_titles';

for i = 1:13
    x{i,3} = [ num2str(round(100*test(1).Model(1).contracts{i}.coinsurance)) '\%' ];
    x{i,5} = [ num2str(round(100*test(1).Model(1).meanCoverage(test(1).Model(1).contracts{i}))) '\%' ];
end
% To print numbers with thousand separator
import java.text.*  
v = DecimalFormat;
for i = 1:13
    x{i,2} = char(v.format(round(test(1).Model(1).contracts{i}.deductible)));
    x{i,4} = char(v.format(round(test(1).Model(1).contracts{i}.oopMax)));
end

fprintf(fid,'\\begin{tabular}{lrrrr}\n');
fprintf(fid,'\\hline\n');
fprintf(fid,'Contract number & Deductible & Coinsurance & Out of pocket maximum & Mean Coverage \\\\ \n');
fprintf(fid,'\\hline\n');  
for i = 1:size(x,1)
    for j = 1:size(x,2)
        if j > 1
           fprintf(fid,' & ');
        end
        fprintf(fid,' %s ',x{i,j});
    end
    fprintf(fid,'\\\\\n');
end
fprintf(fid,'\\end{tabular}');

fclose(fid);

%% Writing results table


fid = fopen('tables/results.tex','w+');
x = cell(13,19);
x(:,1) = contract_titles;
for i = 1:13
    x{i,2} = num2str(round(test(1).DEquilibrium{1}(i),2),'%.2f');
    x{i,8} = num2str(round(test(1).DEfficient{1}(i),2),'%.2f');
    x{i,12} = num2str(round(test(3).DEquilibrium{1}(i),2),'%.2f');
    x{i,18} = num2str(round(test(3).DEfficient{1}(i),2),'%.2f');
    if i > 6
    x{i,5} = num2str(round(test(1).DEquilibrium{2}(i-6),2),'%.2f');
    x{i,15} = num2str(round(test(3).DEquilibrium{2}(i-6),2),'%.2f');
    end
end

% To print numbers with thousand separator
import java.text.*  
v = DecimalFormat;

for i = 1:13
    x{i,3} = char(v.format(round(test(1).pEquilibrium{1}(i))));
    x{i,9} = char(v.format(round(test(1).pEfficient{1}(i))));
    x{i,13} = char(v.format(round(test(3).pEquilibrium{1}(i))));
    x{i,19} = char(v.format(round(test(3).pEfficient{1}(i))));
    if i > 6
    x{i,6} = char(v.format(round(test(1).pEquilibrium{2}(i-6))));
    x{i,16} = char(v.format(round(test(3).pEquilibrium{2}(i-6))));
    end
end

% Printing Headers
fprintf(fid,'\\begin{tabular}{lcclcclccllcclcclcc}\n');
fprintf(fid, '&\\multicolumn{8}{c}{$\\sigma^2_A=0.25$} & & & \\multicolumn{8}{c}{$\\sigma^2_A=2$}\\\\\n');
fprintf(fid, '\\cline{2-9} \\cline{12-19}\n');
fprintf(fid,'& \\multicolumn{2}{c}{No Mandate} &  & \\multicolumn{2}{c}{Mandate} &  & \\multicolumn{2}{c}{Optimum}');
fprintf(fid,'& & & \\multicolumn{2}{c}{No Mandate} &  & \\multicolumn{2}{c}{Mandate} &  & \\multicolumn{2}{c}{Optimum} \\\\\n');
fprintf(fid, '\\cline{2-3} \\cline{5-6} \\cline{8-9} \\cline{12-13} \\cline{15-16} \\cline{18-19}\n');
fprintf(fid, 'Contract number & Q & P & & Q & P & & Q & P & & & Q & P & & Q & P & & Q & P \\\\ \n');
for i = 1:size(x,1)
    for j = 1:size(x,2)
        if j > 1
           fprintf(fid,' & ');
        end
        fprintf(fid,' %s ',x{i,j});
    end
    fprintf(fid,'\\\\\n');
end
fprintf(fid, '&&&&&&&&&&&&&&&&&& \\\\ \n');
fprintf(fid, 'Welfare & \\multicolumn{2}{c}{0} &  & \\multicolumn{2}{c}{%.0f} &  & \\multicolumn{2}{c}{%.0f}',...
    round(test(1).WEquilibrium{2}-test(1).WEquilibrium{1}),round(test(1).WEfficient{1}-test(1).WEquilibrium{1}));
fprintf(fid, '& & & \\multicolumn{2}{c}{0} &  & \\multicolumn{2}{c}{%.0f} &  & \\multicolumn{2}{c}{%.0f}',...
    round(test(3).WEquilibrium{2}-test(3).WEquilibrium{1}),round(test(3).WEfficient{1}-test(3).WEquilibrium{1}));
fprintf(fid,'\\end{tabular}');
fclose(fid);

