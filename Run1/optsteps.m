
function [dist] = optsteps(tx, ty, Sx, Sy, N)
% Calculates optimal path (no diagonals)

dx = abs(tx - Sx); 
dy = abs(ty - Sy); 

d1 = dx + dy; %normal distance
d2 = abs(N - tx + Sx) + dy;  %over the right wall
d3 = dx + abs(N - ty + Sy);  %over top
d4 = abs(N - tx + Sx) + abs(N - ty + Sy); 

DIST = [d1 d2 d3 d4]; 
dist = min(DIST);  

end
