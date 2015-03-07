%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_caralognormal_model_nl_utility %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Eduardo Azevedo and Rafael Mourao
Date:   2015-03-03

This script tests the methods uFunction and cFunction in
the healthcaralognormalmodel_nl class.
%}


clear;
close all;
addpath('../classes')
rng(1);

% Input model parameters
meanS = sqrt(25000^2 - 5100^2);
typeDistributionMean = ...
    [1*10^(-5), 1330, 4340, meanS]; % Original A was 1.9*10^-3
typeDistributionLogCovariance = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
    -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
    -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
    0     0     0    0.25]; % ???

costOfPublicFunds = 0;

% Set up model
deductibleVector  = [1500, 750, 500, 250, 0];
coinsuranceVector = [.1, .1, .1, .1, .1];
oopMaxVector      = [4500, 3750, 3500, 2750, 2500];

Model = healthcaralognormalmodel_nl( ...
    deductibleVector, coinsuranceVector, oopMaxVector, ...
    typeDistributionMean, typeDistributionLogCovariance);

%% Consumer with no moral hazard and small distribution of losses
%{
This agent has no moral hazard, and always loses about $3,000. So he ends
up paying the out of pocket maximum of $1,000, and his gain from insurance
should be $2,000. His cost to the company should be $1,000.

The ex post expenditure should be equal to the loss because there is no
moral hazard.

Ex post utility should go down with a slope of 1 up to the deductible, then
a slope of 0.5, and then be constant once the oop maximum of $1,000 is
reached, which happens at a loss of $1,900.
%}
close all;
x.deductible  = 100;
x.coinsurance = 0.5;
x.oopMax      = 1000;
x.name        = 'test';

type.A = 1e-5;
type.H = 0;
type.M = 3000;
type.S = 10;

display('Contract');
display(x);
display('Agent type: this consumer has a loss of approximately 3,000 and no moral hazard');
display(type);

display('Test exPostUtility');
    nGrid = 30;
    lossGrid = linspace(-500, 3000, nGrid);
    uExPost = zeros(1, nGrid);
    eExPost = zeros(1, nGrid);
    for i = 1:nGrid
        [uExPost(i), eExPost(i), limits] = Model.exPostUtility(x, type, lossGrid(i));
    end;
    display(uExPost);
    display(eExPost);
    
    figure();
    hold on;
    plot(lossGrid, [eExPost; lossGrid]);
    title('Ex post expenditure vs. loss');
    
    figure();
    plot(lossGrid, uExPost);
    title('Ex post utility vs. loss')

display('Test uFunction and cFunction');
    u    = Model.uFunction(x, type);
    c    = Model.cFunction(x, type);
    display(u);
    display(c);
    
%% Consumer with moral hazard and small distribution of losses
%{
This consumer is similar to the previous, still with loss of about $3,000,
but with a moral hazard parameter
of H = $1,000. He will also be in the out of pocket maximum region almost
for sure. He will spend a total of about $4,000, 3 to cover the loss and 1
from moral hazard. Costs to the insurer should be $3,000. The consumer
values this at $2,5000. 3 from the loss, minus 1 from oop, plus 500 from
moral hazard.

ex post expenditure should be $500 for a loss of $0. He always prefers to
go into the copay region, where he spends $500. Expenditures then increase
linearly with the loss (and are always $500 above losses).
Up to some point where he decides to jump over the
oop maximum region (the oop maximum hits at a loss of $1,900). He decides
to jump before hitting $1,900, and jumps to a value that is above $1,900.
from that point on expenditures grow linearly with losses, and are $1,000
above losses (because that is the moral hazard parameter).

ex post utility should be $75 for a loss of $0. He would be happier by
(0.5)^2*H/2= $125 with a linear contract. But in the nonlinear contract he
has to spend an extra $100 (the deductible) minus $50 (how much a linear
contract would charge for that spending). So $125 - $50 = $75. For positive
small losses ex post utility should have a slope of -1. There is then a
threshold where he decides to jump into the out-of-pocket-maximum region.
This jump should coincide with the jump in the expenditure as a function of
losses. In this region ex post utility should be $500 of happiness from the
moral hazard consumption above losses, minus $1,000 from the oop maximum to
a total of -$500.
%}
close all;
type.H = 1000;

display('Contract');
display(x);
display('Agent type: this consumer has a loss of approximately 3,000, but now moral hazard of $1,000');
display(type);

display('Test exPostUtility');
    nGrid = 2000;
    lossGrid = linspace(-500, 3000, nGrid);
    uExPost = zeros(1, nGrid);
    eExPost = zeros(1, nGrid);
    for i = 1:nGrid
        [uExPost(i), eExPost(i), limits] = Model.exPostUtility(x, type, lossGrid(i));
    end;
    display(uExPost);
    display(eExPost);
    
    hold on;
    
    figure();
    plot(lossGrid, [eExPost; lossGrid]);
    title('Ex post expenditure vs. loss')
    
    figure();
    plot(lossGrid, uExPost);
    title('Ex post utility vs. loss')

display('Test uFunction and cFunction');
    u    = Model.uFunction(x, type);
    c    = Model.cFunction(x, type);
    display(u);
    display(c);
    
%% Adding risk aversion.
%{
Now let's make S and A nontrivial, I set it to $10,000 and then $30,000.
To see if we are
calculating things properly with risk. The previous tests suggest we are
doing well ex post, this one is for ex ante.
%}
close all;
type.A = 1.5.*10^(-3);

display('Contract');
display(x);

for i = 1:2
    if i == 1
        type.S = 10000;
    elseif i == 2
        type.S = 30000;
    end;
    
    sprintf('This type has S = %d', type.S)
    display(type);

    display('Test uFunction and cFunction');
        u    = Model.uFunction(x, type);
        c    = Model.cFunction(x, type);
        display(u);
        display(c);
end;

%% Varying risk aversion
%{
This one tests whether the value of a contract goes up sensibly with risk
aversion. Costs should not go up because behavior ex post does not change
with risk aversion.
%}
close all;

display('Contract');
display(x);

nRiskAversionGrid   = 10;
minRiskAversionGrid = 10^(-3);
maxRiskAversionGrid = 10^(-2);
riskAversionGrid    = ...
    linspace(minRiskAversionGrid, maxRiskAversionGrid, nRiskAversionGrid);
uVector = zeros(nRiskAversionGrid, 1);
cVector = zeros(nRiskAversionGrid, 1);

for i = 1:nRiskAversionGrid
    type.A     = riskAversionGrid(i);
    display(type.A);
    uVector(i) = Model.uFunction(x, type);
    cVector(i) = Model.cFunction(x, type);
end;

hold on;

figure();
plot(riskAversionGrid, uVector);
title('Willingness to pay vs risk aversion')

figure();
plot(riskAversionGrid, cVector);
title('Cost vs risk aversion')