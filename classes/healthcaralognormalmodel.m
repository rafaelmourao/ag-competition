classdef healthcaralognormalmodel < model
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
        function Model = healthcaralognormalmodel(slopeVector, ...
                typeDistributionMean, typeDistributionLogCovariance)
            
            Model.typeDistributionMean = typeDistributionMean;
            Model.typeDistributionLogCovariance = ...
                typeDistributionLogCovariance;
            
            n = length(slopeVector);
            for i = 1:n
                x.slope =         slopeVector(i) ;
                x.name  = num2str(slopeVector(i));
                Model.contracts{i} = x;
            end;
        end
            
        function u = uFunction(~, x, type)
                u = x.slope .* type.M ...
                + (x.slope.^2) .* 0.5 .* type.H ...
                + 0.5 .* x.slope .* (2-x.slope) .* (type.S.^2) .* type.A;
        end
        
        function c = cFunction(~, x, type)
                c = x.slope .* type.M ...
                + (x.slope.^2) .* type.H;
        end
        
        function e = eFunction(~, ~, ~)
                e = 0;
        end
            
        function Type = typeDistribution(Model)
            v = Model.lognrndfrommoments(...
                Model.typeDistributionMean, Model.typeDistributionLogCovariance, 1);
            Type.A = v(1);
            Type.H = v(2);
            Type.M = v(3);
            Type.S = v(4);
        end;
        
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

