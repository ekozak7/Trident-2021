function [x] = mod1(number,modnum)
%Modified mod function
%Gets rid of zero

x = mod(number, modnum); 
if x == 0
    x = modnum; 
end

end

