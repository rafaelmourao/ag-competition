classdef model
    %Model Encodes fundamentals of the model.
    %   Instances of this class specify utility and cost functions, a set
    %   of contracts and a distribution of types. Utility and cost
    %   functions are just like in the theoretical model, they take a
    %   contract and type as inputs and produce willingness to pay or
    %   costs. Contracts is a cell with each entry corresponding to a
    %   contract. The user or a subclass can specify these to be anything
    %   as long as the utility and cost functions read them properly.
    %   typeDistribution is a function that returns a randomly drawn type.
    
    properties
        contracts
    end
    
    properties (Dependent)
        nContracts
    end
    
    methods (Abstract)
        u = uFunction(obj,x, type)
        c = cFunction(obj,x, type)
        e = eFunction(obj,x, type) 
        type = typeDistribution(obj)
    end
    
    methods
        % Get Methods
        function n = get.nContracts(obj)
            n = length(obj.contracts);
        end
    end
end

