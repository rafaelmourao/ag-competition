%% Optimize

addpath('../classes');
clear;
close all;
load population.mat;

WEquilibrium = Population.welfare(pEquilibrium, 0);

DEfficient = Population.demand(pEfficient);

f = @(pvar) WEquilibrium - Population.welfare(pvar, 0);
n = 100;

alphagrid = linspace(-0.2, 1, n);
for i = 1:n
    p = alphagrid(i) .* pEfficient + (1-alphagrid(i)) .* pEquilibrium;
    fgrid(i) = f(p);
end;

plot(alphagrid, fgrid);