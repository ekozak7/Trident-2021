function [sx, sy] = movetor(Sx,Sy, direction, N)
%Move in designated direction around a torus
up = [0 1]; 
right = [1 0]; 
down = [0 -1];
left = [-1 0]; 

sx = Sx + direction(1); 
sx = mod(sx, N+1);
if sx == 0
    if direction(1) == 1
        sx = 1; 
    else
        sx = N; 
    end
end

sy = Sy + direction(2); 
sy = mod(sy, N+1);
if sy == 0
    if direction(2) == 1
        sy = 1;
    else
        sy = N; 
    end
end
end

