addpath('../classes')
mkdir('values')
load('tests.mat')

einav_best_contract.deductible = 0;
einav_worst_contract.deductible = 1500;
bronze_contract.deductible = 5200;
einav_best_contract.coinsurance = .1;
einav_worst_contract.coinsurance = .1;
bronze_contract.coinsurance = .35;
einav_best_contract.oopMax = 2500;
einav_worst_contract.oopMax = 3750;
bronze_contract.oopMax = 6400;

actuarial_value_einav_best = round(100*test(1).Model(1).meanCoverage(einav_best_contract));
actuarial_value_einav_worst = round(100*test(1).Model(1).meanCoverage(einav_worst_contract));
actuarial_value_bronze = round(100*test(1).Model(1).meanCoverage(bronze_contract));
actuarial_value_public = round(100*test(1).Model(1).meanCoverage(test(1).Model(1).nullContract));

welfare_gain_optimal_2 = round(test(3).WEfficient{1} - test(3).WEquilibrium{1});
welfare_loss_mandate_2 = round(test(3).WEquilibrium{1} - test(3).WEquilibrium{2});

percentage_bronze_or_less_no_mandate_2 = round(100*sum(test(3).DEquilibrium{1}(1:7)));
percentage_bronze_mandate_2 = round(100*test(3).DEquilibrium{2}(1));

variable_names = {'actuarial_value_einav_best', 'actuarial_value_einav_worst',...
    'actuarial_value_bronze', 'actuarial_value_public', 'welfare_gain_optimal_2',...
    'welfare_loss_mandate_2', 'percentage_bronze_or_less_no_mandate_2',...
    'percentage_bronze_mandate_2'};

for i = variable_names
    fid = fopen(['values/' i{:} '.tex'],'w+');
    fprintf(fid,'%.0f',eval(i{:}));
    fclose(fid);
end