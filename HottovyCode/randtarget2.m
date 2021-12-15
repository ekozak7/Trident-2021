function [tx, ty] = randtarget2(N)
%Places a target randomly within a NxN grid
%Prevents zero index (target along axes)

tx = randi(N);
ty = randi(N); 

end

