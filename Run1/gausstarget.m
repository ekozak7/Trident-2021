function [tx, ty] = gausstarget(N, sigma)
%Places a target according to the Gaussian distribution within a NxN grid
%Prevents zero index (target along axes) 

m = N/2; 

tx = round(normrnd(m, sigma)); 
ty = round(normrnd(m, sigma)); 

tx = mod(tx, N);
if tx == 0
    tx = N; 
end
ty = mod(ty, N); 
if ty == 0
    ty = N; 
end

end