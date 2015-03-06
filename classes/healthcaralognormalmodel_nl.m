classdef healthcaralognormalmodel_nl < model
    %   healthcaralognormalmodel_nl Creates a health insurance model as in
    %   the Azevedo and Gottlieb paper. Subclass of the model class.
    %   This model has CARA preferences, normal losses, lognormally
    %   distributed types, and contracts are parametrized by three elements
    %   as in Einav et al. (2013). Inputs to the constructor:
    %   deductibleVector is a vector of deductibles, coinsuranceVector is a vector
    %   of coinsurance shares, oopVector is a vector with out of pocket maximum
    %   values, parameterMean is a 4 dimensional vector of means of parameters
    %   and parameterLogVariance is a 4x4 matrix of log covariances.
    %   Parameters are ordered as A, H, M, S as in the Azevedo and Gottlieb
    %   paper (absolute risk aversion A, moral hazard H, mean loss M,
    %   and standard deviation of losses S).
    
    properties
        nullContract
        typeDistributionMean
        typeDistributionLogCovariance
    end
    
    methods
        % Constructor
        function obj = healthcaralognormalmodel_nl(deductibleVector, ...
                coinsuranceVector, oopMaxVector, typeDistributionMean, typeDistributionLogCovariance)
            
            obj.typeDistributionMean = typeDistributionMean;
            obj.typeDistributionLogCovariance = typeDistributionLogCovariance;
            n = length(deductibleVector);
            for i = 1:n
                x.deductible       = deductibleVector(i) ;
                x.coinsurance      = coinsuranceVector(i) ;
                x.oopMax           = oopMaxVector(i) ;
                x.name             = num2str(i);
                obj.contracts{i} = x;
            end;
            obj.nullContract.deductible  = Inf;
            obj.nullContract.coinsurance = 1;
            obj.nullContract.oopMax      = Inf;
            obj.nullContract.name        = 'Null Contract';
        end
        
        function u = uFunction(obj, x, type)
            
            [u0, ~, bounds] = exPostUtility(obj, x, type, 0);
            limits = integrationLimits(obj, type, 1e-2);
            u = 0;
            
            if (limits(1) < 0)
                
                u = integral(@(x) lossDistributionFunction(obj, type, x), ...
                    limits(1), 0,'AbsTol', 1e-15,'RelTol',1e-12 )...
                    * -exp(-type.A * u0);
                
            end
                
            
            u = u + integral(@(l) -exp(-type.A*exPostUtility(obj, x, type, l)) .*...
                lossDistributionFunction(obj, type, l), 0, limits(2),...
                'AbsTol', 1e-15,'RelTol',1e-12,'WayPoints',bounds(isfinite(bounds)));
            
            % Eduardo's addition. Calculate utility from no insurance.
            
            u0 = integral(@(l) -exp(-type.A * ...
                exPostUtility(obj, obj.nullContract, type, l)) .*...
                lossDistributionFunction(obj, type, l),max(0,limits(1)),limits(2),...
                'AbsTol', 1e-15,'RelTol',1e-12);
            
            if (u0 - u > 1e-6)
                error('Utility without insurance cannot be higher than with it')
            end

            % Calculate certainty equivalent
            
            CE  = log(u ./ u0) ./ (-type.A);
            u   = CE;
            
        end
        
        function c = cFunction(obj, x, type)
            
            [c0,~] = exPostCost(obj, x, type, 0);
            limits = integrationLimits(obj, type, 1e-2);
            c = 0;
            
            if (limits(1) < 0)
                
                c = integral(@(x) lossDistributionFunction(obj,type,x), ...
                    limits(1),max(0,limits(2)),...
                    'AbsTol', 1e-15,'RelTol',1e-12)*c0;
                
            end
            
            if (limits(2) > 0)
            
            c = c + integral(@(l) exPostCost(obj, x, type, l).*...
                lossDistributionFunction(obj, type, l),max(0,limits(1)),limits(2),...
                'AbsTol', 1e-15,'RelTol',1e-12);
            
            end
            
        end
        
        function Type = typeDistribution(obj)
            
            v = lognrndfrommoments(...
                obj.typeDistributionMean, obj.typeDistributionLogCovariance, 1);
            
            Type.A = v(1);
            Type.H = v(2);
            Type.M = v(3);
            Type.S = v(4);
            
            function v = ...
                    lognrndfrommoments(meanVector, logCovMatrix, varargin)
                %   estimate_lognormal_from_moments Estimate lognormal parameters
                %   Inputs: a line vector of means, and a log covariance matrix.
                %   Outputs: A line vector of draws. Optional argument for number of lines.
                %   Warning: for some reason MATLAB chokes if one of the variables has 0
                %   variance, so always put at least a tiny bit of bogus variance.
                
                nDraws = 1;
                if (length(varargin) == 1)
                    nDraws = varargin{1};
                end;
                
                % Find lognormal parameters
                sigmaMatrix = cholcov(logCovMatrix);
                b  = sum(sigmaMatrix.^2, 1);
                mu = log(meanVector) - b/2;
                
                % Draw
                x = randn(nDraws, length(meanVector));
                v = exp(repmat(mu,[nDraws, 1]) + x * sigmaMatrix);
            end
            
        end
        
        function limits = integrationLimits(obj, type, tol)
            % Finds the smallest possible interval such that all points in
            % it have positive (at machine accuracy) density. Perfoms a
            % check to see if integrating the density between these points
            % results in 1
            
            limits(1) = findCloserZero(type.M,1e6,tol);
            limits(2) = findCloserZero(type.M,-1e6,tol);
            
            integralCheck = integral(@(l)...
                lossDistributionFunction(obj, type, l), ...
                limits(1),limits(2),'AbsTol', 1e-15,'RelTol',1e-12);
            
            if ( abs(integralCheck - 1) > 1e-6 )
                error('Integral could not be well approximated, or this is not a distribution')
            end
            
            function b = findCloserZero(b,d_init,tol)
                f_b = 0;
                d = d_init;
                while (abs(d) > tol || f_b > 0)
                    b = b - d;
                    f_b = obj.lossDistributionFunction(type,b);
                    if (f_b == 0)
                        if (sign(d) == sign(d_init))
                            d = - d / 10;
                        end
                    else
                        if (sign(d) ~= sign(d_init))
                            d = -d / 10;
                        end
                    end
                end
            end
            
        end
        
        function x = lossDistributionFunction(~,type,l)
            x = normpdf(l,type.M,type.S);
        end
        
        function [cost, bounds] = exPostCost(obj, x, type, losses)
            [~, e, c, bounds] = exPostUtility(obj, x, type, losses);
            cost = e - c;            
        end
        
        function [u, e, c, bounds] = exPostUtility(~, x, type, losses)
            u = zeros(1, length(losses));
            e = zeros(1, length(losses));
            c = zeros(1, length(losses));
            bounds = zeros(length(losses), 3);
            for i = 1:length(losses)
                l = max(losses(i), 0);
                bounds(i, 1) = max(min(x.deductible-(1-x.coinsurance)*type.H/2,x.oopMax-type.H/2),0);
                bounds(i, 2) = max(x.deductible-(1-x.coinsurance)*type.H/2,0);
                bounds(i, 3) = max((x.oopMax-(1-x.coinsurance)*x.deductible)/x.coinsurance ...
                    - (2 - x.coinsurance) * type.H / 2, 0);
                if l < bounds(i, 1)
                    u(i)= -l;
                    e(i) = l;
                    c(i) = l;
                elseif (l >= bounds(i, 2)) && (l <= bounds(i, 3))
                    u(i) = (1-x.coinsurance)^2*type.H/2 - (1-x.coinsurance)*x.deductible - x.coinsurance*l;
                    e(i) = (1-x.coinsurance)*type.H + l;
                    c(i) = x.deductible + (1-x.coinsurance)*(e(i)-x.deductible);
                else
                    u(i) = type.H/2 - x.oopMax;
                    e(i) = type.H + l;
                    c(i) = x.oopMax;
                end
            end
        end
        
        function [populationSize, CalculationParametersEquilibrium, CalculationParametersOptimum] = ...
                suggestComputationParameters(Model, percentError)
            
            nContracts = Model.nContracts;
            Population = population(Model, 100);
            priceOrderOfMagnitude = mean(Population.cMatrix(:));
            
            populationSize = ...
                floor(nContracts / percentError^2 / 2);
            CalculationParametersEquilibrium.behavioralAgents = ...
                percentError;
            CalculationParametersEquilibrium.fudge = ...
                percentError / nContracts / 100;
            CalculationParametersEquilibrium.tolerance = ...
                percentError * priceOrderOfMagnitude;
            CalculationParametersEquilibrium.maxIterations = ...
                floor(10 / CalculationParametersEquilibrium.fudge);
            
            CalculationParametersOptimum.tolerance = CalculationParametersEquilibrium.tolerance;
            CalculationParametersOptimum.maxIterations = 10^4;
        end;
        
    end
end

