function [sx,sy] = randmovetor(Sx,Sy, N)
%randomly move on a torus

y = randi(2);
y = mod(y, 2);
yp = rand - .5;
yp = yp/abs(yp);
y = y * yp;

x = y - 1; 
x = mod(x, 2); 
xp = rand - .5; 
xp = xp/abs(xp); 
x = x*xp; 

direction = [x, y]; 

sx = Sx + direction(1); 
sx = mod(sx + N+1, N+1); 

sy = Sy + direction(2); 
sy = mod(sy + N+1, N+1);
end

