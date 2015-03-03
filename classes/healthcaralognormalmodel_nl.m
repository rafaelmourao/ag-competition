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
        typeDistributionMean
        typeDistributionLogCovariance
    end
    
    methods
        % Constructor
        function Model = healthcaralognormalmodel_nl(deductibleVector, ...
                coinsuranceVector, oopMaxVector, typeDistributionMean, typeDistributionLogCovariance)
            
            Model.typeDistributionMean = typeDistributionMean;
            Model.typeDistributionLogCovariance = typeDistributionLogCovariance;
            n = length(deductibleVector);
            for i = 1:n
                x.deductible       = deductibleVector(i) ;
                x.coinsurance      = coinsuranceVector(i) ;
                x.oopMax           = oopMaxVector(i) ;
                x.name             = num2str(i);
                Model.contracts{i} = x;
            end;
        end
        
        function u = uFunction(obj, x, type)
            
            [u, ~, limits] = exPostUtility(obj, x, type, 0);
            u = integral(@(x) lossDistributionFunction(obj, type, x), -Inf, 0) * u;
            if limits(1) > 0
                u = u + integral(@(l) -exp(-type.A * exPostUtility(obj, x, type, l)).*...
                    lossDistributionFunction(obj, type, l), 0, limits(1));
            end
            if limits(3) > limits(2)
                u = u + integral(@(l) -exp(-type.A*exPostUtility(obj, x, type, l)) .*...
                    lossDistributionFunction(obj,type,l),limits(2),limits(3));
            else
                u = u + integral(@(l) -exp(-type.A*exPostUtility(obj, x, type, l)) .*...
                    lossDistributionFunction(obj,type,l),limits(1),Inf);
                return
            end
            if ~isinf(limits(3))
                u = u + integral(@(l) -exp(-type.A * exPostUtility(obj, x, type, l)) .*...
                    lossDistributionFunction(obj, type, l),limits(3),Inf);
            end
            
            % Eduardo's addition. Calculate utility from no insurance.
            xNull.deductible  = Inf;
            xNull.coinsurance = 1;
            xNull.oopMax      = Inf;
            xNull.name        = 'Null Contract';
            u0 = integral( ...
                @(l) lossDistributionFunction(obj, type, l) ...
                .* -exp(-type.A * exPostUtility(obj, x, type, l)) ...
                , -Inf, 0);
            
            % Calculate certainty equivalent
            CE  = log(u ./ u0) ./ (-type.A);
            u   = CE;
        end
        
        function u = uFunction_alt(obj, x, type)
            
            [u, ~, limits] = exPostUtility(obj, x, type, 0);
            u = integral(@(x) lossDistributionFunction(obj, type, x), -Inf, 0) * u;
            u = u + integral(@(l) -exp(-type.A*exPostUtility(obj, x, type, l)) .*...
                    lossDistributionFunction(obj, type, l), 0, Inf);
        end
        
         function u = uFunction_alt2(obj, x, type)
            
            [u, ~, limits] = exPostUtility(obj, x, type, 0);
            u = integral(@(x) lossDistributionFunction(obj, type, x), -Inf, 0) * u;
            if limits(1) > 0
                u = u + integral(@(l) -exp(-type.A*exPostUtility_alt(obj, x, type, l)).*...
                    lossDistributionFunction(obj,type,l),0,limits(1));
            end
            if limits(3) > limits(2)
                u = u + integral(@(l) -exp(-type.A*exPostUtility_alt(obj, x, type, l)).*...
                    lossDistributionFunction(obj,type,l),limits(2),limits(3));
            else
                u = u + integral(@(l) -exp(-type.A*exPostUtility_alt(obj, x, type, l)).*...
                    lossDistributionFunction(obj,type,l),limits(1),Inf);
                return
            end
            if ~isinf(limits(3))
                u = u + integral(@(l) -exp(-type.A*exPostUtility_alt(obj, x, type, l)).*...
                    lossDistributionFunction(obj,type,l),limits(3),Inf);
            end     
            
         end
        
         function u = uFunction_alt3(obj, x, type)
            
            [u, ~, limits] = exPostUtility(obj, x, type, 0);
            u = integral(@(x) lossDistributionFunction(obj,type,x),-Inf,0)*u;
            u = u + integral(@(l) -exp(-type.A*exPostUtility_alt(obj, x, type, l)).*...
                    lossDistributionFunction(obj,type,l),0,Inf);
        end
        
        function c = cFunction(obj, x, type)
            
            [~,c,limits] = exPostUtility(obj, x, type, 0);
            c = integral(@(x) lossDistributionFunction(obj,type,x),-Inf,0)*c;
            if limits(1) > 0
                c = c + integral(@(l) exPostExpenditure(obj, x, type, l).*...
                    lossDistributionFunction(obj,type,l),0,limits(1));
            end
            if limits(3) > limits(2)
                c = c + integral(@(l) exPostExpenditure(obj, x, type, l).*...
                    lossDistributionFunction(obj,type,l),limits(2),limits(3));
            else
                c = c + integral(@(l) exPostExpenditure(obj, x, type, l).*...
                    lossDistributionFunction(obj,type,l),limits(1),1e6);
                return
            end
            if ~isinf(limits(3))
                c = c + integral(@(l) exPostExpenditure(obj, x, type, l).*...
                    lossDistributionFunction(obj,type,l),limits(3),1e6);
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
        
        function x = lossDistributionFunction(~,type,l)
            x = normpdf(l,type.M,type.S);
        end
        
        function [e, limits] = exPostExpenditure(obj, x, type, losses)
            [~, e, limits] = exPostUtility(obj, x, type, losses);
        end
        
        function [u, e, limits] = exPostUtility(~, x, type, losses)
            u = zeros(1, length(losses));
            e = zeros(1, length(losses));
            limits = zeros(length(losses), 3);
            for i = 1:length(losses)
                l = max(losses(i), 0);
                limits(i, 1) = max(min(x.deductible-(1-x.coinsurance)*type.H/2,x.oopMax-type.H/2),0);
                limits(i, 2) = max(min(x.deductible-(1-x.coinsurance)*type.H/2,0));
                limits(i, 3) = max((x.oopMax-(1-x.coinsurance)*x.deductible)/x.coinsurance ...
                    - (2 - x.coinsurance) * type.H / 2, 0);
                if l < limits(i, 1)
                    u(i)= -l;
                    e(i) = l;
                elseif (l >= limits(i, 2)) && (l <= limits(i, 3))
                    u(i) = (1-x.coinsurance)^2*type.H/2 - (1-x.coinsurance)*x.deductible - x.coinsurance*l;
                    e(i) = (1-x.coinsurance)*type.H + l;
                else
                    u(i) = type.H/2 - x.oopMax;
                    e(i) = type.H + l;
                end
            end
        end
        
        function [u, e] = exPostUtility_alt(~, x, type, l)
            l = max(l,0);
            u(1,:) = -l;
            e(1,:) = l;
            u(2,:) = (1-x.coinsurance)^2*type.H/2 - (1-x.coinsurance)*x.deductible - x.coinsurance*l;
            e(2,:) = (1-x.coinsurance)*type.H + l;
            u(3,:) = type.H/2 - x.oopMax;
            e(3,:) = type.H + l;
            [u, loc] = max(u);
            e = e(loc);
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

