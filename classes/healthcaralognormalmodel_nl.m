classdef healthcaralognormalmodel_nl < model
    %   healthcaralognormalmodel_nl Creates a health insurance model as in
    %   the Azevedo and Gottlieb paper. Subclass of the model class.
    %   This model has CARA preferences, lognormal losses, lognormally
    %   distributed types, and contracts are parametrized by three elements
    %   as in Einav et al. (2013). Inputs to the constructor:
    %   deductibleVector is a vector of deductibles, coinsuranceVector is a
    %   vector of coinsurance shares, oopMaxVector is a vector with out of
    %   pocket maximum values, publicInsuranceMaximum is the maximum of
    %   losses an uninsured individual must cover before public insurance
    %   public insurance covers the rest, parameterMean is a 4 dimensional
    %   vector of means of parameters and parameterLogVariance is a 4x4
    %   matrix of log covariances.
    %   Parameters are ordered as A, H, M, S as in the Azevedo and Gottlieb
    %   paper (absolute risk aversion A, moral hazard H, mean loss M,
    %   and standard deviation of losses S).
    
    properties
        publicInsuranceMaximum
        nullContract
        typeDistributionMean
        typeDistributionLogCovariance
    end
    
    methods ( Access = public, Hidden = false ) % Public methods
        
        function obj = healthcaralognormalmodel_nl(deductibleVector, ...
                coinsuranceVector, oopMaxVector, publicInsuranceMaximum, ...
                typeDistributionMean, typeDistributionLogCovariance)
            
            obj.typeDistributionMean = typeDistributionMean;
            obj.typeDistributionLogCovariance = typeDistributionLogCovariance;
            n = length(deductibleVector);
            for i = 1:n
                contract.deductible       = deductibleVector(i) ;
                contract.coinsurance      = coinsuranceVector(i) ;
                contract.oopMax           = oopMaxVector(i) ;
                contract.name             = num2str(i);
                obj.contracts{i} = contract;
            end;
            
            obj.publicInsuranceMaximum   = publicInsuranceMaximum;
            obj.nullContract.deductible  = publicInsuranceMaximum;
            obj.nullContract.coinsurance = 1;
            obj.nullContract.oopMax      = publicInsuranceMaximum;
            obj.nullContract.name        = 'Null Contract';
        end
        
        function x = lossPDF(~, type, l)
            % Probability distribution function for the losses
            sigma = sqrt( log( 1 + (type.S./type.M).^2 ) );
            mu    = log(type.M) - sigma.^2 / 2;
            x = lognpdf(l, mu, sigma);
        end
        
        function x = lossCDF(~, type, l)
            % Cumulative distribution function for the losses
            sigma = sqrt( log( 1 + (type.S./type.M).^2 ) );
            mu    = log(type.M) - sigma.^2 / 2;
            x = logncdf(l, mu, sigma);
        end
        
        function u = uFunction(obj, contract, type)
            % Returns the willingness to pay of an agent for a contract,
            % considering the outside option of public insurance.
            
            % Checking contract
            checkContract(obj, contract);
            
            % If contract is null, then return zero
            if (contract.deductible == obj.publicInsuranceMaximum)
                u = 0;
                return
            end
            
            if (type.A > 0)
                
                u = logExpectedExponentialValue( obj, @(l) -type.A*...
                    exPostUtility(obj, contract, type, l),  contract, type );
                
                u0 = logExpectedExponentialValue( obj, @(l) -type.A*...
                    exPostUtility(obj, obj.nullContract, type, l), ...
                    obj.nullContract, type );
                
                u = - u / type.A;
                u0 = - u0 / type.A;
                
            else % If there is no risk aversion, just calculate the 
                 %  expected ex post utility
                
                u = expectedValue( obj, @(l) ...
                    exPostUtility(obj, contract, type, l),  contract, type );
                
                u0 = expectedValue( obj, @(l) ...
                    exPostUtility(obj, obj.nullContract, type, l), ...
                    obj.nullContract, type );
                
            end
            
            if ( (u0 - u)/abs(u) > 1e-6 )
                fprintf('Utility with insurance: %.8f\n',u)
                fprintf('Utility without insurance: %.8f\n',u0)
                error('Utility without insurance cannot be higher than with it')
            end
            
            u  = u - u0;
            
            if ( u > obj.publicInsuranceMaximum + type.H / 2 )
                u = obj.publicInsuranceMaximum + type.H / 2;
                warning('Impossible value, setting u = publicInsuranceMaximum + type.H / 2')
            end
            
        end
        
        function c = cFunction(obj, contract, type)
            % Returns the expected private cost incurred by the insurer
            % when the agent chooses a contract. Currently, costs
            % occur only when the agent chooses private insurance
            
            % Checking contract
            checkContract(obj, contract);
            
            % If contract is null, then return zero cost
            if (contract.deductible == obj.publicInsuranceMaximum)
                c = 0;
            else
                c = expectedValue( obj, @(l) ...
                    exPostCost(obj, contract, type, l), contract, type );
            end
            
        end
        
        function e = eFunction(obj, contract, type)
            % Returns the expected externality generated when the
            % agent chooses a contract. Currently, externalities occur
            % only when the agent chooses public insurance.
            
            % Checking contract
            checkContract(obj, contract);
            
            if (contract.deductible ~= obj.publicInsuranceMaximum)
                e = 0;
            else
                e = expectedValue( obj, @(l) ...
                    (l>obj.publicInsuranceMaximum).*(l-obj.publicInsuranceMaximum),...
                    contract, type );
            end
            
        end
        
        function type = typeDistribution(obj)
            % Returns a rando structure with the individual caracteristics
            % distributed according to the lognormal parameters
            
            v = lognrndfrommoments(...
                obj.typeDistributionMean, obj.typeDistributionLogCovariance, 1);
            
            type.A = v(1);
            type.H = v(2);
            type.M = v(3);
            type.S = v(4);
            
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
        
        function m = meanCoverage(obj, contract)
            % Returns the expected coverage for the mean individual
            
            checkContract(obj, contract);
            
            meantype.A = obj.typeDistributionMean(1);
            meantype.H = obj.typeDistributionMean(2);
            meantype.M = obj.typeDistributionMean(3);
            meantype.S = obj.typeDistributionMean(4);
            
            mexpenditure =  expectedValue( obj, @(l) ...
                exPostExpenditure(obj, contract, meantype, l), contract, meantype );
            
            if (contract.deductible == obj.publicInsuranceMaximum)
                
                mcost = eFunction(obj, contract, meantype);
                
            else
                
                mcost = cFunction(obj, contract, meantype);
                
            end
            
            m = mcost/mexpenditure;
            
        end
        
        function expenditure = exPostExpenditure(obj, contract, type, losses)
            
            if (contract.deductible == obj.publicInsuranceMaximum)
                expenditure = losses;
            else
                [~, expenditure] = exPostUtility(obj, contract, type, losses);
            end
            
        end
        
        function cost = exPostCost(obj, contract, type, losses)
            
            if (contract.deductible == obj.publicInsuranceMaximum)
                cost = 0;
            else
                [~, expenditure, payment] = exPostUtility(obj, contract, type, losses);
                cost = expenditure - payment;
            end
            
        end
        
        function [u, expenditure, payment, bounds] = exPostUtility(obj, contract, type, losses)
            
            % Initializing variables
            u = zeros(1, length(losses));
            expenditure = zeros(1, length(losses));
            payment = zeros(1, length(losses));
            bounds = zeros(length(losses), 3);
            
            for i = 1:length(losses)
                
                % Loss is bounded below by zero
                l = max(losses(i), 0);
                
                % If contract is the null contract
                if (contract.deductible == obj.publicInsuranceMaximum)
                    
                    bounds(i,:) = obj.publicInsuranceMaximum;
                    u(i) = max(-l,-obj.publicInsuranceMaximum);
                    payment(i) = min(l,obj.publicInsuranceMaximum);
                    expenditure(i) = l;
                    
                else
                    
                    % Calculating the loss boundaries for each interval of expenses
                    bounds(i, 1) = max(min(contract.deductible-(1-contract.coinsurance)*type.H/2,contract.oopMax-type.H/2),0);
                    bounds(i, 2) = max(contract.deductible-(1-contract.coinsurance)*type.H/2,0);
                    bounds(i, 3) = max((contract.oopMax-(1-contract.coinsurance)*contract.deductible)/contract.coinsurance ...
                        - (2 - contract.coinsurance) * type.H / 2, 0);
                    
                    if l < bounds(i, 1)
                        
                        expenditure(i) = l;
                        u(i) = -l;
                        payment(i) = l;
                        
                    elseif (l >= bounds(i, 2)) && (l < bounds(i, 3))
                        
                        u(i) = (1-contract.coinsurance)^2*type.H/2 - (1-contract.coinsurance)*contract.deductible - contract.coinsurance*l;
                        expenditure(i) = (1-contract.coinsurance)*type.H + l;
                        payment(i) = contract.deductible + contract.coinsurance*(expenditure(i)-contract.deductible);
                        
                    else
                        
                        u(i) = type.H/2 - contract.oopMax;
                        expenditure(i) = type.H + l;
                        payment(i) = contract.oopMax;
                        
                    end
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
        end
    end
    
    methods ( Access = private, Hidden = true ) % Auxiliary methods
        
        function checkContract(obj, contract)
            if ( contract.coinsurance < 0 || contract.coinsurance > 1 )
                error('Coinsurance must be between zero and one')
            elseif ( contract.deductible > contract.oopMax || contract.deductible < 0 )
                error(['Deductible must be higher than zero, and lower than'...
                    ' the out of pocket maximum'])
            elseif ( contract.oopMax > obj.publicInsuranceMaximum )
                error(['Out of pocket maximum must be lower than the public'...
                    'insurance maximum'])
            end
        end
        
        function x = logExpectedExponentialValue(obj, function_handle, contract, type )
            
            [limits, oopMaxLoss] = integralLimits(obj, contract, type);
            
            if ( limits(1) ~= limits(2) )
                [~, K] = fminbnd(@(l) - log( lossPDF(obj, type, l) )...
                    - function_handle(l), limits(1), limits(2) );
                K = -K;
            else
                K = 0;
            end
            
            x = leftIntegral(obj,  @(l) exp( function_handle(l) - K ), ...
                contract, type, limits );
            
            x = x + innerIntegral (obj, @(l) exp( log( lossPDF(obj, type, l) ) ...
                + function_handle(l) - K), contract, type, limits, oopMaxLoss );
            
            x = x + rightIntegral( obj,@(l) exp( function_handle(l) - K ),...
                contract, type, limits, oopMaxLoss );
            
            x = log(x) + K;
            
        end
        
        function x = expectedValue(obj, function_handle, contract, type )
            
            
            [limits, oopMaxLoss] = integralLimits(obj, contract, type);
            
            x = leftIntegral(obj, function_handle, contract, type, limits );
            
            x = x + innerIntegral (obj, @(l) lossPDF(obj, type, l)...
                .* function_handle(l), contract, type, limits, oopMaxLoss );
            
            x = x + rightIntegral(obj, @(l) function_handle(l), contract,...
                type, limits, oopMaxLoss );
            
        end
        
        function x = leftIntegral(obj, function_handle, ~, type, limits )
            
            x = 0;
            
            if ( limits(1) == 0 )
                x = lossCDF(obj, type, 0) * function_handle(0);
            end
            
        end
        
        function x = innerIntegral(obj, function_handle, contract, type, limits, oopMaxLoss )
            
            x = 0;
            
            [~, ~, ~, bounds] = exPostUtility(obj, contract, type, 0);
            
            if ( limits(2) > 0 || limits(1) <  oopMaxLoss )
                x  = integral(@(l) function_handle(l), limits(1), limits(2),...
                    'AbsTol', 1e-15,'RelTol',1e-12,'WayPoints',...
                    [bounds(isfinite(bounds)),linspace(limits(1), limits(2),1e3)] );
            end
            
        end
        
        function x = rightIntegral(obj, function_handle, contract, type, limits, oopMaxLoss )
            
            x = 0;
            
            if ( limits(2) == oopMaxLoss )
                
                inclination = function_handle( oopMaxLoss + 1 )...
                    - function_handle( oopMaxLoss );
                
                x = (1 - lossCDF(obj, type, oopMaxLoss)) ...
                    * function_handle( oopMaxLoss );
                
                if (inclination > 0)
                    
                    x = x + inclination...
                        * ( type.M - innerIntegral(obj,@(l) l.*lossPDF(obj, type, l),...
                        contract, type, limits, oopMaxLoss) -  (1 - lossCDF(obj, type, oopMaxLoss)) ...
                        * oopMaxLoss  );
                    
                end
                
            end
            
        end
        
        function [limits, oopMaxLoss] = integralLimits(obj, contract, type)
            
            [~, ~, ~, bounds] = exPostUtility(obj, contract, type, 0);
            oopMaxLoss = max(bounds(1),bounds(3));
            
            if (lossPDF(obj,type,0) > 0);
                limits(1) = 0;
            elseif ( (type.M > oopMaxLoss ) && ...
                    lossPDF(obj, type, oopMaxLoss) == 0)
                limits(1) = oopMaxLoss;
            else
                limits(1) = max(0,findCloserNonZero(min(type.M,oopMaxLoss),type.M/100,1e-10));
            end
            
            if (lossPDF(obj,type,oopMaxLoss) > 0);
                limits(2) = oopMaxLoss;
            elseif (type.M > oopMaxLoss )
                limits(2) = oopMaxLoss;
            else
                limits(2) = max(0, findCloserNonZero(type.M,-type.M/100,1e-10));
            end
            
            function b = findCloserNonZero(b,d_init,tol)
                f_b = 0;
                d = d_init;
                while ( abs(d) > tol || f_b == 0 )
                    b = b - d;
                    f_b = obj.lossPDF(type,b);
                    if (f_b > 0)
                        if (sign(d) ~= sign(d_init))
                            d = - d / 10;
                        end
                    else
                        if (sign(d) == sign(d_init))
                            d = - d / 10;
                        end
                    end
                end
            end
        end
    end
    
end

