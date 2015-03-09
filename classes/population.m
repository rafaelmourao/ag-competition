classdef population
    %population A finite list of simulated consumers.
    %   This class encodes a finite list of consumer types and a utility and cost
    %   matrices associated with this population and a model.
    %   These matrices specify willingness to pay (cost) of each consumer i for product j.
    %   The
    %   constructor takes the model and the population size as inputs. This
    %   class also contains methods for calculating demand, equilibrium,
    %   and optimal allocations.
    
    properties
        typeList
        uMatrix
        cMatrix
        size
        nContracts
    end
    
    methods
        % Contructor
        % This constructor is pretty slow right now, because it loops over
        % each agent in the population. I did it like this to make it
        % easier to use different subclasses of the model class. But it is
        % super under-vectorized, so could be optimized a lot at the cost
        % of some flexibility.
        function Population = population(Model, n, nworkers)
            if ( nargin < 3 || (nworkers < 0) )
                nworkers = 0;
            end
            n_Contracts = Model.nContracts;
            u_Matrix = zeros(n, Model.nContracts);
            c_Matrix = zeros(n, Model.nContracts);
            type_List = cell(1,n);
            parfor(i = 1 : n, nworkers)
                type = typeDistribution(Model);
                type_List{i} = type;
                for j = 1 : n_Contracts
                    x = Model.contracts{j};
                    u_Matrix(i, j) = Model.uFunction(x, type);
                    c_Matrix(i, j) = Model.cFunction(x, type);
                end
            end
            Population.uMatrix    = u_Matrix;
            Population.cMatrix    = c_Matrix;
            Population.typeList   = type_List;
            Population.size       = n;
            Population.nContracts = n_Contracts;
        end
        
        % Basic functions
        
        function [D, TC, CS, choiceVector] = demand(Population, p)
            % demand: Takes as input the population and a price vector.
            % Outputs demand vector, total cost vector, and consumer suplus
            % vectors. These are 1xnContracts sized vectors. Also outputs a
            % choiceVector which is populationSize x 1, and each entry
            % specifies the contract j chosen by consumer i.
            
            % surplus = Population.uMatrix - repmat(p, Population.size, 1);
            % Alternative way, faster for large populations:
            surplus = bsxfun(@minus,Population.uMatrix,p);
            [maxSurplus, choiceVector] = max(surplus, [], 2);
            
            % Getting the matrix indices for the maximum surpluses
            % choiceIndices = sub2ind(size(surplus),(1:Population.size)',choiceVector)
            % alternative, faster way:
            choiceIndices =  (1:Population.size)' + (choiceVector-1)*Population.size;
            
            % Using histc to get the frequency of every contract choice
            D = histc(choiceVector',1:Population.nContracts)/Population.size;
            
            % Calculating the cost for each contract by using the
            % matlab function accumarray, summing all expected costs 
            % according to the choice vector subscripts
            TC = accumarray(choiceVector, Population.cMatrix(choiceIndices), ...
                [Population.nContracts,1])'/Population.size;
            
            if nargout > 2
            CS = mean(maxSurplus);
            end
        end
                       
        function W = welfare(Population, p, costOfPublicFunds)
            [D, TC, CS, ~] = Population.demand(p);
            W = CS + (1+costOfPublicFunds).*(D * p' - sum(TC));
        end;      
        
        % Computational methods
        function [p, D, AC, ComputationOutput] = findequilibrium(Population, CalculationParameters)
            % findequilibrium: Finds an equilibrium by iterating average cost. To make it
            % numerically stable must move towards average cost only a
            % small fraction of the way. Parameters:
            % CalculationParameters.behavioralAgents,
            % CalculationParameters.fudge (what fraction of the way you move towards average cost,
            % 1 being the fastest and something close to 0 being more numerically stable,
            % CalculationParameters.maxIterations, CalculationParameters.tolerance
            % Output: price, demand, average cost and a ComputationOutput
            % struct with fields .nIterations, .error (mean square distance
            % between p and AC. Note that this number is often big even
            % with a precise computation because it is driven by contracts
            % taht are not traded. A good improvement would be having a
            % better definition of error in each step of the algorithm),
            % .runTime.
            % This numerical method is pretty accurate for finding
            % equilibrium prices within the range of contracts that are
            % actually traded in equilibrium. The prices of contracts that
            % are not traded are less stable. They become even less stable
            % in a model with a lot of contracts. Note that the fudge
            % factor has to be pretty low to compute things accurately
            % because of the issue of contracts that are not traded. When a
            % contract just starts being traded the AC curve can be very
            % steep, and the algorithm will not converge. The behavioral
            % agents attenuate this tendency, so a small number of
            % behavioral agents make the numerics more finicky. The reason
            % for using this method as opposed to something smarter that
            % changes fudge factors etc based on the AC curves is that this
            % is simple to implement, and typically works if we set the
            % computational parameters conservatively. But there is a lot
            % of room for improvement.
            tic;
            epsilon = CalculationParameters.behavioralAgents / Population.nContracts; % epsilon is the number of behavioral agents per contracts.
            fudge   = CalculationParameters.fudge;
            
            error       = Inf;
            nIterations = 0;
            
            function [p1, error] = iteration(p0)
                [D, TC] = Population.demand(p0);
                AC            = TC ./ (D + epsilon);
                error         = norm(AC - p0, Inf);
                currentFudge  = fudge + 1.1^(-nIterations-1); % I made the fudge factor close to 1 in the first few iterations so that it moves fast in the beginning. But decreasing by 10% in each iteration so that it quickly gets to the value specified in the function call.
                p1            = currentFudge*AC + (1-currentFudge)*p0;
            end
            
            p = zeros(1, Population.nContracts);
            while (error > CalculationParameters.tolerance) ...
                    && (nIterations < CalculationParameters.maxIterations)
                [p, error] = iteration(p);
                % Require at least 50 iterations.
                if (nIterations < 50)
                    error = Inf;
                end;
                nIterations = nIterations + 1;
            end;
            
            ComputationOutput.error       = error;
            ComputationOutput.nIterations = nIterations;
            ComputationOutput.runTime     = toc;
        end
                
        function [p, W, ComputationOutput] = findefficient(Population, costOfPublicFunds, CalculationParameters)
            tic;
            % findefficient: This function finds an efficient allocation
            % given a cost of public funds. Inputs a population, cost of
            % public funds, and calculation parameters. This is a
            % structures with fields .tolerance and .maxIterations. Outputs
            % are prices, welfare, and a ComputationOutput struct with
            % fields .nIterations, .error (in % terms vis a vis last
            % iteration), and .runTime. The algorithm is designed to work
            % with ordered sets of contracts. Each iteration loops over all
            % contracts. It tries to maximize the price difference between
            % contract j and j+1 while keeping fixed the relative prices
            % of contracts below j and above j+1. So the numerical method
            % is similar to the perturbation proof approach in Saez (2002).
            
            % Calculate maximum and minimum marginal utility.
            MU = diff(Population.uMatrix, 1, 2);
            
            % Use this to create upper and lower bound for derivatives and initial condition.
            dp_max = max(MU);
            lower_bound = dp_max .* 0;
            upper_bound = dp_max .* 1.2;
            dp0 = (lower_bound+upper_bound)/2;
            
            % Define function that turns dp into p, function that evaluates
            % -welfare given dp (the difference vector of p), and function that updates the ith
            % coordinate if dp (so that we can choose each dp(j) at a time to maximize welfare).
            integratedp = @(dp) cumsum([0, dp]);
            f = @(dp)  -Population.welfare(integratedp(dp), costOfPublicFunds);
            function dpOut = dp_update(dpjIn, jIn, dpIn)
                dpOut      = dpIn;
                dpOut(jIn) = dpjIn;
            end
            
            
            % Initialize loop variables
            nIterations = 0;
            error = Inf;
            W0 = Inf;
            dp = dp0;
            
            % Main loop.
            while (nIterations < CalculationParameters.maxIterations) ...
                    && (error > CalculationParameters.tolerance)
                % Update each dp(j) maximizing as in the Saez (2002) perturbation.
                for j = 1 : length(dp)
                    g = @(dpj) f( dp_update(dpj,j,dp) );
                    [dpj, W] = fminbnd(g, lower_bound(j), upper_bound(j));
                    dp(j) = dpj;
                end;
                
                % Update iteration, error, and improvement in welfare.
                nIterations = nIterations + 1;
                if nIterations > 50 % Guarantee that at least 50 iterations are done.
                    error = norm(W-W0) ./ norm(W);
                end;
                W0 = W;
            end;
            
            p = integratedp(dp);
            W = -W;
            
            ComputationOutput.nIterations = nIterations;
            ComputationOutput.error       = error;
            ComputationOutput.runTime     = toc;
        end;
        
        % Graphing methods
        function graphHandle = graphEFC(Population)
            % graphEFC: This function plots graphs similar to the Einav
            % Finkelstein and Cullen paper. These are for models with two contracts.
            % Plots of average cost difference between contracts, the
            % demand for the high contract as a function of the price
            % difference, and the marginal total cost of increasing
            % coverage. These graphs are similar to EFC when one of the
            % contracts corresponds to buying nothing. They get a little
            % weird when one of the contracts corresponds to some positive
            % level of coverage. See the Veiga and Weyl Leibniz rule paper
            % for details on this. Inputs are a population, and the output
            % is a figure handle.
            if Population.nContracts ~= 2
                graphHandle = -99;
                display('Must use model with two contracts.');
            else
                nPointstoPlot = 100;
                wtp = diff(Population.uMatrix, 1, 2);
                dpVector = linspace(min(wtp), max(wtp), nPointstoPlot);
                for i = 1 : nPointstoPlot-1
                    dp = dpVector(i);
                    p = [0, dp];
                    [D, TC, ~, ~] = Population.demand(p);
                    pVector(i)   = p(2);
                    qVector(i)   = D(2);
                    dACVector(i) = TC(2)/D(2) - TC(1)/D(1);
                    TCVector(i)  = sum(TC);
                end;
                
                I = (qVector > 0.10) & (qVector<1);
                
                pFit     = fit(qVector', pVector'  , 'smoothingspline');
                ACFit    = fit(qVector(I)', dACVector(I)', 'smoothingspline');
                TCFit    = fit(qVector', TCVector' , 'poly3');
                
                %                 ACdiff = differentiate(ACFit, qVector);
                %                 MCVector = qVector .* ACdiff' + feval(ACFit, qVector)';
                TCdiff = differentiate(TCFit, qVector);
                MCVector =TCdiff';
                
                graphHandle = figure();
                hold on;
                plot(qVector(I), pVector(I)           );
                plot(qVector(I), dACVector(I), 'red'  );
                plot(qVector(I), MCVector(I) , 'green');
            end;
        end;
    end
    
end

