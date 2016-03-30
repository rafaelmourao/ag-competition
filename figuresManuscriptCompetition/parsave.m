function parsave(fname, varnames, varargin)
numvars=numel(varargin);
fname = [fname, '.mat'];

if (length(varnames) ~= numvars)
    display(['Error! Vector of names has a different number of ' ...
             'elements than the number of variables!'])

for i=1:numvars
       eval([varnames(i), '=varargin{i};']);  
end

save('-mat', fname, varnames(1));

for i = 2:numvars    
    save('-mat', fname, varnames(i), '-append');
end
end