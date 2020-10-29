function [x] = drawLevy_Hottovy(c,mu)
%Draw a sample from the Levy distribution
%parameters c and mu

    x = random('stable',mu-1,1,c,c,1);

end

