%MCTS known, multiple targets
%on a grid, torus
%non destructive search

clear all, clc

%% Parameters

T = 1;  %Number of trials
N = 20;    %Size of search grid
M = N + 1; %Size of matrices
L = 100;  %Number of loops
c = 2;   %exploration constant
rv = 1; %vision radius
Limit = 10000; 
Steps = 100; 

targs = 2; 
Tx = zeros(1, targs); 
Ty = zeros(T, targs); 

%Vectors
    up = [0 1]; 
    right = [1 0]; 
    down = [0 -1];
    left = [-1 0];


%% Data Collection

TG = zeros(1, T); 

%% Run Trials
tic

for t = 1:T
   
     t
     j = 0; %step count
     
    %Store Data
    Nx = zeros(M); %visit count
    Qx = zeros(M); %success rate
    K = ones(M);  %total step count
    Kx = zeros(M); %average step count
    
    %Place Targets
    for i = 1:targs
        Tx(i) = randi(N); 
        Ty(i) = randi(N); 
    end
    
    %Place Searcher
    Sx = N/2; 
    Sy = N/2;
    
    %Record Path Traveled
    Px = [Sx]; 
    Py = [Sy]; 
    
    %Count Targets
    Targets = 0; 
    X = [];
    Y = []; 
    
    %Start MCTS
    
    D = zeros(1, targs);
    for d = 1:targs
        D(d) = dist(Tx(d), Ty(d), Sx, Sy); 
    end
    [mindist, place] = min(D); 
    
    while j <= Steps
       
        if Nx(Sx+1, mod1(Sy + 2,M)) == 0  
            [sx, sy] = movetor(Sx, Sy, up, N);  %move up
            Nx(Sx+1, Sy+1) = Nx(Sx+1, Sy+1) + 1;  %add visit count to root node
            Nx(Sx+1, mod1(Sy+2,M)) = Nx(Sx+1, mod1(Sy+2,M)) + 1; %add visit count to up node
            sx1 = sx; %create disposable variables
            sy1 = sy; 
            count =0; 
            D = zeros(1, targs);
               for d = 1:targs
                   D(d) = dist(Tx(d), Ty(d), sx1, sy1); 
               end
               [mindist, place] = min(D); 
             
             while mindist > rv %random playout
                    count = count + 1; 
                    if count > Limit
                        fprintf('limit! 1\n')
                        break
                    end
                [sx1, sy1] = randmovetor(sx1, sy1, N);
                K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %add step count to selected node
                Nx(sx1+1, sy1+1) = Nx(sx1+1, sy1+1) + 1; %add visit count to playout nodes
               %Check min distance
                D = zeros(1, targs);
               for d = 1:targs
                   D(d) = dist(Tx(d), Ty(d), sx1, sy1); 
               end
               [mindist, place] = min(D);
             end
            Kx(sx+1, sy+1) = K(sx+1, sy+1)/Nx(sx+1, sy+1);  %calculate avg step count
            %Qx(sx+1, sy+1) = 1/Kx(sx+1, sy+1); %update win count
            if Kx(sx+1, sy+1) + 1 >= max(Kx)
                Qx(sx+1, sy+1) = (1+ Targets)/Kx(sx+1, sy+1);
            else
            Qx(sx+1, sy+1) = (1)/Kx(sx+1, sy+1);
            end
       end
        %Second node: Right
        if Nx(mod1(Sx+2,M), Sy + 1) == 0  
            [sx, sy] = movetor(Sx, Sy, right, N);  %move right
            Nx(Sx+1, Sy+1) = Nx(Sx+1, Sy+1) + 1;  %add visit count to root node
            Nx(mod1(Sx+2,M), Sy+1) = Nx(mod1(Sx+2, M), Sy+1) + 1; %add visit count to up node
            sx1 = sx; %create disposable variables
            sy1 = sy; 
            count = 0; 
                D = zeros(1, targs);
               for d = 1:targs
                   D(d) = dist(Tx(d), Ty(d), sx1, sy1); 
               end
               [mindist, place] = min(D); 
             
             while mindist > rv %random playout
                    count = count + 1; 
                    if count > Limit
                        fprintf('limit! 1\n')
                        break
                    end
                [sx1, sy1] = randmovetor(sx1, sy1, N); 
                K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %add step count to selected node
                Nx(sx1+1, sy1+1) = Nx(sx1+1, sy1+1) + 1; %add visit count to playout nodes
                 D = zeros(1, targs);
               for d = 1:targs
                   D(d) = dist(Tx(d), Ty(d), sx1, sy1); 
               end
               [mindist, place] = min(D); 
                end
            Kx(sx+1, sy+1) = K(sx+1, sy+1)/Nx(sx+1, sy+1);  %calculate avg step count
            %Qx(sx+1, sy+1) = 1/Kx(sx+1, sy+1);     %update win count
           if Kx(sx+1, sy+1) + 1 >= max(Kx)
                Qx(sx+1, sy+1) = (1+ Targets)/Kx(sx+1, sy+1);
            else
            Qx(sx+1, sy+1) = (1)/Kx(sx+1, sy+1);
            end
        end
        %Third node: Down
        if Nx(Sx+1, mod1(Sy ,M)) == 0  
            [sx, sy] = movetor(Sx, Sy, down, N);  %move down
            Nx(Sx+1, Sy+1) = Nx(Sx+1, Sy+1) + 1;  %add visit count to root node
            Nx(Sx+1, mod1(Sy,M)) = Nx(Sx+1, mod1(Sy,M)) + 1; %add visit count to up node
            sx1 = sx; %create disposable variables
            sy1 = sy; 
            count = 0; 
                D = zeros(1, targs);
               for d = 1:targs
                   D(d) = dist(Tx(d), Ty(d), sx1, sy1); 
               end
               [mindist, place] = min(D); 
             
             while mindist > rv %random playout
                    count = count + 1; 
                    if count > Limit
                        fprintf('limit! 1\n')
                        break
                    end
                [sx1, sy1] = randmovetor(sx1, sy1, N); 
                K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %add step count to selected node
                Nx(sx1+1, sy1+1) = Nx(sx1+1, sy1+1) + 1; %add visit count to playout nodes
             D = zeros(1, targs);
               for d = 1:targs
                   D(d) = dist(Tx(d), Ty(d), sx1, sy1); 
               end
               [mindist, place] = min(D);    
             end
            Kx(sx+1, sy+1) = K(sx+1, sy+1)/Nx(sx+1, sy+1);  %calculate avg step count
            %Qx(sx+1, sy+1) = 1/Kx(sx+1, sy+1);     %update win count
            if Kx(sx+1, sy+1) + 1 >= max(Kx)
                Qx(sx+1, sy+1) = (1+ Targets)/Kx(sx+1, sy+1);
            else
            Qx(sx+1, sy+1) = (1)/Kx(sx+1, sy+1);
            end
        end
        %Fourth node: Left
        if Nx(mod1(Sx, M), Sy + 1) == 0  
            [sx, sy] = movetor(Sx, Sy, left, N);  %move left
            Nx(Sx+1, Sy+1) = Nx(Sx+1, Sy+1) + 1;  %add visit count to root node
            Nx(Sx, Sy+1) = Nx(Sx, Sy+1) + 1; %add visit count to up node
            sx1 = sx; %create disposable variables
            sy1 = sy; 
            count = 0; 
               D = zeros(1, targs);
               for d = 1:targs
                   D(d) = dist(Tx(d), Ty(d), sx1, sy1); 
               end
               [mindist, place] = min(D); 
             
             while mindist > rv %random playout
                    count = count + 1; 
                    if count > Limit
                        fprintf('limit! 1\n')
                        break
                    end
                [sx1, sy1] = randmovetor(sx1, sy1, N); 
                K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %add step count to selected node
                Nx(sx1+1, sy1+1) = Nx(sx1+1, sy1+1) + 1; %add visit count to playout nodes
             D = zeros(1, targs);
               for d = 1:targs
                   D(d) = dist(Tx(d), Ty(d), sx1, sy1); 
               end
               [mindist, place] = min(D);    
             end
            Kx(sx+1, sy+1) = K(sx+1, sy+1)/Nx(sx+1, sy+1);  %calculate avg step count
            %Qx(sx+1, sy+1) = 1/Kx(sx+1, sy+1);     %update win count
            if Kx(sx+1, sy+1) + 1 >= max(Kx)
                Qx(sx+1, sy+1) = (1+ Targets)/Kx(sx+1, sy+1);
            else
            Qx(sx+1, sy+1) = (1)/Kx(sx+1, sy+1);
            end
        end
      
       %Run simulations
       for loops = 1:L
           
          sx = Sx;  %reset starting node 
          sy = Sy;
          Nx(sx+1, sy+1) = Nx(sx+1, sy+1) + 1;  %add visit count to root node
          UCT = Qx./Nx + c*sqrt(abs(log(Nx(Sx+1, Sy+1))./Nx)); %calculate UCT
          %restrict to four movements
          UCTi = zeros(1, 4); 
              UCTi(1) = UCT(sx+1, mod1(sy+2,M)); %up
              UCTi(2) = UCT(mod1(sx+2, M), sy+1); %right 
              UCTi(3) = UCT(sx+1, mod1(sy, M)); %down
              UCTi(4) = UCT(mod1(sx, M), sy+1); %left
          [Max, R] = max(UCTi);  %Max is max, R is direction
          %move in max direction
          if R == 1
              [sx,sy] = movetor(sx,sy, up, N); 
          elseif R == 2
              [sx,sy] = movetor(sx,sy, right, N);
          elseif R == 3
              [sx,sy] = movetor(sx,sy, down, N);
          else
              [sx,sy] = movetor(sx,sy, left, N);
          end
          %disposable variables
          sx1 = sx; 
          sy1 = sy; 
          Nx(sx+1, sy+1) = Nx(sx+1, sy+1) + 1;  %add visit count for selected node
          %random playout
          count = 0; 
          D = zeros(1, targs);
               for d = 1:targs
                   D(d) = dist(Tx(d), Ty(d), sx1, sy1); 
               end
               [mindist, place] = min(D); 
             
             while mindist > rv %random playout
                    count = count + 1; 
                    if count > Limit
                        fprintf('limit! 1\n')
                        break
                    end
              [sx1, sy1] = randmovetor(sx1, sy1, N); 
              K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %add step count
              Nx(sx1+1, sy1+1) = Nx(sx1+1, sy1+1) + 1; %add visit count to playout nodes
          D = zeros(1, targs);
               for d = 1:targs
                   D(d) = dist(Tx(d), Ty(d), sx1, sy1); 
               end
               [mindist, place] = min(D); 
             end
          Kx(sx+1, sy+1) = K(sx+1, sy+1)/Nx(sx+1, sy+1); %calculate avg step count
          %Qx(sx+1, sy+1) = 1/Kx(sx+1, sy+1);  %calculate win count
          if Kx(sx+1, sy+1) + 1 >= max(Kx)
                Qx(sx+1, sy+1) = (1+ Targets)/Kx(sx+1, sy+1);
            else
            Qx(sx+1, sy+1) = (1)/Kx(sx+1, sy+1);
            end
       end
       
       %Play best move
            UCT = Qx./Nx + c*sqrt(abs(log(Nx(Sx+1, Sy+1))./Nx)); %calculate UCT
            %restrict to four moves
            UCTi = zeros(1, 4);
                UCTi(1) = UCT(Sx+1, mod1(Sy+2, M)); %up
                UCTi(2) = UCT(mod1(Sx+2, M), Sy+1); %right 
                UCTi(3) = UCT(Sx+1, mod1(Sy, M)); %down
                UCTi(4) = UCT(mod1(Sx, M), Sy+1); %left
                
                UCTi
          [Max, R] = max(UCTi);  %Max is max, R is direction
          %move in max direction
          if R == 1
              [Sx,Sy] = movetor(Sx,Sy, up, N); 
          elseif R == 2
              [Sx,Sy] = movetor(Sx,Sy, right, N);
          elseif R == 3
              [Sx,Sy] = movetor(Sx,Sy, down, N);
          else
              [Sx,Sy] = movetor(Sx,Sy, left, N);
          end
          %Add to path
          Px = [Px Sx];
          Py = [Py Sy] ;
          
          j = length(Px); 
          
          %Check for targets
          D = zeros(1, targs);
               for d = 1:targs
                   D(d) = dist(Tx(d), Ty(d), Sx, Sy); 
                   if D(d) <= rv
                       Targets = Targets +1 
                       X = [X Tx(d)];
                       Y = [Y Ty(d)]; 
                   end
               end
               
      
    end
    
    sametarg = 0; 
   for j = 2:Targets
       if (X(j) == X(j-1)) && (Y(j) == Y(j-1))
           sametarg = sametarg + 1; 
       end
   end
   
   Targets = Targets - sametarg; 
   
   TG(t) = Targets; 
     
    
end

time = toc

Avgtarg = mean(TG)
Vartarg = var(TG)
Avgtime = time/T

figure
for i = 1:targs
   plot(Tx(i), Ty(i), '*r')
   hold on  
end
plot(Px, Py, '-b')
axis([0 N 0 N])
tit = sprintf('%d Targets found after %d Steps', Targets, length(Px))
title(tit);
