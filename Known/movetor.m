function [sx, sy] = movetor(Sx,Sy, direction, N)
%Move in designated direction around a torus

up = [0 1]; 
right = [1 0]; 
down = [0 -1];
left = [-1 0];

if direction == 1
    direction = [0 1]; 
elseif direction == 2
    direction = [1 0]; 
elseif direction == 3
    direction = [0 -1]; 
else
    direction = [-1 0]; 
end
 

sx = Sx + direction(1); 
sx = mod(sx, N+1);
if sx == 0
    if direction(1) == -1
        sx = N; 
    else
        sx = 1; 
    end
end

sy = Sy + direction(2); 
sy = mod(sy, N+1);
if sy == 0
    if direction(2) == -1
        sy = N; 
    else 
        sy = 1; 
    end
end
end

