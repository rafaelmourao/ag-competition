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
        uFunction
        cFunction
        typeDistribution
    end
    
    properties (Dependent)
        nContracts
    end
    
    methods
        % Empty constructor
        function Model = model()
        end;
        
        % Get Methods
        function n = get.nContracts(Model)
            n = length(Model.contracts);
        end;
    end
end

