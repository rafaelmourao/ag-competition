function [str]=num2bank(num)
% num2bank      Turns number into financial string
%       Donwloaded off internet and slightly modified. My version does not
%       display cents and puts a dollar sign in front.
     str = arrayfun(@(x) num2bankScalar(x) , num, 'UniformOutput', false); 
end

function [str]=num2bankScalar(num)
     num=floor(num*100)/100;
     str = num2str(num);
     k=find(str == '.', 1);
     if(isempty(k))
         str=[str,'.00'];
     end
     FIN = min(length(str),find(str == '.')-1);
     for i = FIN-2:-3:2
         str(i+1:end+1) = str(i:end);
         str(i) = ',';
     end
     
     str = ['$', str];
     str = str(1:length(str)-3);
end