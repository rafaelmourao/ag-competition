classdef healthcaralognormalmodel_nl < model
    %healthcaralognormalmodel Creates a health insurance model as in
    %the Azevedo and Gottlieb paper. Subclass of the model class.
    %   This models have CARA preferences, normal losses, linear contracts,
    %   and lognormally distributed types. Inputs to the constructor:
    %   slopeVector is a vector of contract slopes, parameterMean is a 4
    %   dimensional vector of means of parameters and parameterLogVariance
    %   is a 4x4 matrix of log covariances. Parameters are ordered as A, H,
    %   M, S as in the Azevedo and Gottlieb paper (absolute risk aversion A
    %   , moral hazard H, mean loss M, and standard deviation of losses S).
    
    properties
        typeDistributionMean
        typeDistributionLogCovariance
    end
    
    methods
        % Constructor
        function Model = healthcaralognormalmodel_nl(deductibleVector, ...
                coinsuranceVector, oopVector, typeDistributionMean, typeDistributionLogCovariance)
            Model.typeDistributionMean = typeDistributionMean;
            Model.typeDistributionLogCovariance = typeDistributionLogCovariance;
            
            n = length(deductibleVector);
            for i = 1:n
                x.deductible =         deductibleVector(i) ;
                x.coinsurance =        coinsuranceVector(i) ;
                x.oop =                oopVector(i) ;
                x.name  = num2str(i);
                Model.contracts{i} = x;
            end;
        end
            
        function u = uFunction(obj, x, type)
            [u,~,limits] = exPostUtility(obj, x, type, 0);
            u = integral(@(x) lossDistributionFunction(obj,type,x),-Inf,0)*u;
            if limits(1) > 0 
                u = u + integral(@(l) -exp(-type.A*exPostUtility(obj, x, type, l))*...
                    lossDistributionFunction(obj,type,l),0,limits(1));
            end
            if limits(3) > limits(2)
                u = u + integral(@(l) -exp(-type.A*exPostUtility(obj, x, type, l))*...
                    lossDistributionFunction(obj,type,l),limits(2),limits(3));
            else
                u = u + integral(@(l) -exp(-type.A*exPostUtility(obj, x, type, l))*...
                    lossDistributionFunction(obj,type,l),limits(1),Inf);
                return
            end
            u = u + integral(@(l) -exp(-type.A*exPostUtility(obj, x, type, l))*...
                    lossDistributionFunction(obj,type,l),limits(3),Inf);
        end
        
        function c = cFunction(~, x, type)
            [c,limits] = exPostUtility(obj, x, type, 0);
            c = integral(@(x) lossDistributionFunction(obj,type,x),-Inf,0)*c;
            if limits(1) > 0 
                c = c + integral(@(l) exPostExpenditure(obj, x, type, l)*...
                    lossDistributionFunction(obj,type,l),0,limits(1));
            end
            if limits(3) > limits(2)
                c = c + integral(@(l) exPostExpenditure(obj, x, type, l)*...
                    lossDistributionFunction(obj,type,l),limits(2),limits(3));
            else
                c = c + integral(@(l) exPostExpenditure(obj, x, type, l)*...
                    lossDistributionFunction(obj,type,l),limits(1),Inf);
                return
            end
                c = c + integral(@(l) exPostExpenditure(obj, x, type, l)*...
                    lossDistributionFunction(obj,type,l),limits(3),Inf);
        end
            
        function Type = typeDistribution(Model)
            v = Model.lognrndfrommoments(...
                Model.typeDistributionMean, Model.typeDistributionLogCovariance, 1);
            Type.A = v(1);
            Type.H = v(2);
            Type.M = v(3);
            Type.S = v(4);
        end;
        
        function [e, limits] = exPostExpenditure(obj, x, type, l)
            [~,e,limits] = exPostUtility(obj, x, type, l);
        end
        
        function [u, e, limits] = exPostUtility(~, x, type, l)
            l = max(l,0);
            limits(1) = max(min(x.deductible-(1-x.coinsurance)*type.H/2,x.oop-type.H/2),0);
            limits(2) = max(min(x.deductible-(1-x.coinsurance)*type.H/2,0));
            limits(3) = max((x.oop-(1-x.coinsurance)*x.deductible)/x.coinsurance ...
                    - (1-x.coinsurance)*type.H/2,0);
            if l < limits(1)
                u = -l;
                e = l;
            elseif (l >= limits(2)) && (l <= limits(3))
                u = (1-x.coinsurance)^2*type.H/2 - (1-x.coinsurance)*x.deductible - x.coinsurance*l;
                e = (1-x.coinsurance)*type.H + l;
            else
                u = type.H/2 - x.oop;
                e = type.H + l;
            end
        end
        
        function [u, e] = exPostUtility_alt(~, x, type, l)
            l = max(l,0);
            u(1) = -l;
            e(1) = l;
            u(2) = (1-x.coinsurance)^2*type.H/2 - (1-x.coinsurance)*x.deductible - x.coinsurance*l;
            e(2) = (1-x.coinsurance)*type.H + l;
            u(3) = type.H/2 - x.oop;
            e(3) = type.H + l;
            [u, loc] = max(u);
            e = e(loc);
        end
         
        function x = lossDistributionFunction(~,type,l)
            x = normcdf(type.M,type.S,l);
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
    
    methods(Access = private, Static = true)       
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
        end;
    end
end

