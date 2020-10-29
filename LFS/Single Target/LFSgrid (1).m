%Levy on a grid, torus

clear all, clc

T = 1000;  %Number of trials
N = 40;    %Size of search grid
rv = 1; 
c = 1; 
mu = 2 ; 

Steps = zeros(1, T); 
Optimal = zeros(1, T); 
Diff = zeros(1, T); 

up = [0 1]; 
right = [1 0]; 
down = [0 -1];
left = [-1 0];

%sigma = 10;

tic 

for t = 1:T
    t
   %[Tx, Ty] = gausstarget(N, sigma); 
   Tx = randi(N); 
   Ty = randi(N); 
   
   Sx = 1; 
   Sy = 1; 
   
   Px = [];
   Py = []; 
   
   Optimal(t) = optsteps(Tx, Ty, Sx, Sy, N); 
   
   while dist(Tx, Ty, Sx, Sy) > rv
       
       theta = rand*2*pi;  
       h = drawLevy_Hottovy(c,mu);
       
       % Check if flight is too long or too short
            count = 0;
            while log10(h)>4.5 || h<rv
                h = drawLevy_Hottovy(c,mu);
                count = count+1;
                %sprintf('Too long of flight!!!')
                if count > 100
                    sprintf('Whoa!!!!')
                    count
                end
            end
            
            x = round(h*cos(theta)); 
            y = round(h*sin(theta)); 
           
            sx1 = mod(Sx + x, N+1);  %line will be wrong if this is greater than N
            sy1 = mod(Sy + y, N+1);
            
            %equation for line
            %m = (sy1 - Sy)/(sx1 - Sx); 
            %lX = linspace(Sx, sx1); 
            %lY = m*lX - m*Sx + Sy; 
            
            %closest point
                %P = [Sx, (Sy + 1); (Sx + 1), Sy; Sx, (Sy -1); (Sx - 1), Sy];
                %IO = [sx1, sy1]; 
                %u = [(sx1 - Sx), (sy1 - Sy)]; 
                %[d2H, H] = point_to_line_distance(P, u, IO)
               count = 0;  
            while (Sx ~= sx1) || (Sy ~= sy1)
               U = dist(sx1, sy1, Sx, Sy+1); 
               R = dist(sx1, sy1, Sx+1, Sy); 
               D = dist(sx1, sy1, Sx, Sy-1);
               L = dist(sx1, sy1, Sx-1, Sy); 
               O = [U, R, D, L]; 
               [M, A] = min(O);
               O = [up; right; down; left];
               B = O(A,:);
               [Sx, Sy] = movetor(Sx, Sy, B, N);
               Px = [Px Sx];
               Py = [Py Sy];
               count = count + 1;
               if dist(Sx, Sy, Tx, Ty) <= rv
                   break
               end
              
                
            end  
   end
   
   Steps(t) = length(Px);
   Diff(t) = Steps(t) - Optimal(t); 
    
end

step = Steps(t); 
AvgSteps = mean(Steps)
Variance = var(Steps)

time = toc

AvgDiff = mean(Diff) 
VarDiff = var(Diff)

%figure
%plot(Px, Py, '-b')
%hold on
%plot(Tx, Ty, '*r')
%axis([0 20 0 20])
%t = sprintf('%d Steps', step);
%title(t)
%print('lfsg8', '-dpng');


