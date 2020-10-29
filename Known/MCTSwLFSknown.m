%MCTS with LFS Direction bias
%one Target
%Searcher starts in corner
%Target is known
%Torus

%uses LFS for distance, MCTS for direction

clear all, clc

%% Parameters

T = 1000;  %Number of trials
N = 20;    %Size of search grid
M = N+1;
L = 10;  %Number of loops
c = 2;   %exploration constant
rv = 1; 
mu = 2; 
Limit = 10000; 

up = [0 1]; 
right = [1 0]; 
down = [0 -1];
left = [-1 0];

sigma = 100; 

%% Data Collection

MCTS = zeros(1, T);  %Number of steps per trial
MOPT = zeros(1, T);  %Optimal steps per trial

MAvg = mean(MCTS);   %Average steps per trial
MVar = var(MCTS);    %Variance of steps per trial

Mdiff = MCTS - MOPT; %Difference from optimal
MAdiff = mean(Mdiff); %Average distance from optimal
MVardiff = var(Mdiff); %Variance of difference from optimal

%% Run Trials

tic
for t = 1:T
   t
    %Store Data
    Nx = zeros(M); %visit count
    Qx = zeros(M); %success rate
    K = zeros(M);  %total step count
    Kx = zeros(M); %average step count
    
    %Vectors
    up = [0 1]; 
    right = [1 0]; 
    down = [0 -1];
    left = [-1 0];

    %Place Target
    [Tx, Ty] = randtarget2(N);  
    
    %Place Searcher
    Sx = 1; 
    Sy = 1;
    
    %Record Path Traveled
    Px = [Sx]; 
    Py = [Sy]; 
    
    %Measure Optimal Path
    MOPT(t) = optsteps(Tx, Ty, Sx, Sy, N); 
    
    %Start MCTS
    while dist(Sx, Sy, Tx, Ty) > rv
    %while (Sx~= Tx) || (Sy ~= Ty)
        if length(Px) > Limit
            Px;
            break
            fprintf('Error'); %output something
        end
       
        %%Test unvisited nodes
        %First node: UP
       if Nx(Sx+1, mod1(Sy + 2,M)) == 0  
            [sx, sy] = movetor(Sx, Sy, up, N);  %move up
            Nx(Sx+1, Sy+1) = Nx(Sx+1, Sy+1) + 1;  %add visit count to root node
            Nx(Sx+1, mod1(Sy+2,M)) = Nx(Sx+1, mod1(Sy+2,M)) + 1; %add visit count to up node
            sx1 = sx; %create disposable variables
            sy1 = sy; 
            count =0; 
            
            h = drawLevy_Hottovy(c,mu);
            h = round(h); 
            dh = 0; 
                while dist(Tx, Ty, sx1, sy1) > rv%random playout
                    count = count + 1; 
                    if count > Limit
                        fprintf('limit! 1\n')
                        break
                    end
                    if dh < h
                        [sx1, sy1] = movetor(sx1, sy1, up, N); 
                        dh = dh + 1; 
                    else
                    [sx1, sy1] = randmovetor(sx1, sy1, N);
                    end
                K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %add step count to selected node
                Nx(sx1+1, sy1+1) = Nx(sx1+1, sy1+1) + 1; %add visit count to playout nodes
                end
            Kx(sx+1, sy+1) = K(sx+1, sy+1)/Nx(sx+1, sy+1);  %calculate avg step count
            Qx(sx+1, sy+1) = 1/Kx(sx+1, sy+1);     %update win count
       end
        %Second node: Right
        if Nx(mod1(Sx+2,M), Sy + 1) == 0  
            [sx, sy] = movetor(Sx, Sy, right, N);  %move right
            Nx(Sx+1, Sy+1) = Nx(Sx+1, Sy+1) + 1;  %add visit count to root node
            Nx(mod1(Sx+2,M), Sy+1) = Nx(mod1(Sx+2, M), Sy+1) + 1; %add visit count to up node
            sx1 = sx; %create disposable variables
            sy1 = sy; 
            count = 0; 
            h = drawLevy_Hottovy(c,mu);
            h = round(h); 
            dh = 0; 
                while dist(Tx, Ty, sx1, sy1) > rv%random playout
                    count = count + 1; 
                    if count > Limit
                        fprintf('limit! 1\n')
                        break
                    end
                    if dh < h
                        [sx1, sy1] = movetor(sx1, sy1, right, N); 
                        dh = dh + 1; 
                    else
                    [sx1, sy1] = randmovetor(sx1, sy1, N);
                    end
                K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %add step count to selected node
                Nx(sx1+1, sy1+1) = Nx(sx1+1, sy1+1) + 1; %add visit count to playout nodes
                end
            Kx(sx+1, sy+1) = K(sx+1, sy+1)/Nx(sx+1, sy+1);  %calculate avg step count
            Qx(sx+1, sy+1) = 1/Kx(sx+1, sy+1);     %update win count
        end
        %Third node: Down
        if Nx(Sx+1, mod1(Sy ,M)) == 0  
            [sx, sy] = movetor(Sx, Sy, down, N);  %move down
            Nx(Sx+1, Sy+1) = Nx(Sx+1, Sy+1) + 1;  %add visit count to root node
            Nx(Sx+1, mod1(Sy,M)) = Nx(Sx+1, mod1(Sy,M)) + 1; %add visit count to up node
            sx1 = sx; %create disposable variables
            sy1 = sy; 
            count = 0;
            h = drawLevy_Hottovy(c,mu);
            h = round(h); 
            dh = 0; 
                while dist(Tx, Ty, sx1, sy1) > rv%random playout
                    count = count + 1; 
                    if count > Limit
                        fprintf('limit! 1\n')
                        break
                    end
                    if dh < h
                        [sx1, sy1] = movetor(sx1, sy1, down, N); 
                        dh = dh + 1; 
                    else
                    [sx1, sy1] = randmovetor(sx1, sy1, N);
                    end
                [sx1, sy1] = randmovetor(sx1, sy1, N); 
                K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %add step count to selected node
                Nx(sx1+1, sy1+1) = Nx(sx1+1, sy1+1) + 1; %add visit count to playout nodes
                end
            Kx(sx+1, sy+1) = K(sx+1, sy+1)/Nx(sx+1, sy+1);  %calculate avg step count
            Qx(sx+1, sy+1) = 1/Kx(sx+1, sy+1);     %update win count
        end
        %Fourth node: Left
        if Nx(mod1(Sx, M), Sy + 1) == 0  
            [sx, sy] = movetor(Sx, Sy, left, N);  %move left
            Nx(Sx+1, Sy+1) = Nx(Sx+1, Sy+1) + 1;  %add visit count to root node
            Nx(Sx, Sy+1) = Nx(Sx, Sy+1) + 1; %add visit count to up node
            sx1 = sx; %create disposable variables
            sy1 = sy; 
            count = 0;
            h = drawLevy_Hottovy(c,mu);
            h = round(h); 
            dh = 0; 
                while dist(Tx, Ty, sx1, sy1) > rv%random playout
                    count = count + 1; 
                    if count > Limit
                        fprintf('limit! 1\n')
                        break
                    end
                    if dh < h
                        [sx1, sy1] = movetor(sx1, sy1, left, N); 
                        dh = dh + 1; 
                    else
                    [sx1, sy1] = randmovetor(sx1, sy1, N);
                    end 
                K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %add step count to selected node
                Nx(sx1+1, sy1+1) = Nx(sx1+1, sy1+1) + 1; %add visit count to playout nodes
                end
            Kx(sx+1, sy+1) = K(sx+1, sy+1)/Nx(sx+1, sy+1);  %calculate avg step count
            Qx(sx+1, sy+1) = 1/Kx(sx+1, sy+1);     %update win count
        end
      
       %Run simulations
       for loops = 1:L
           
          sx = Sx;  %reset starting node 
          sy = Sy;
          Nx(sx+1, sy+1) = Nx(sx+1, sy+1) + 1;  %add visit count to root node
          
          
          %% Choose direction based on UCT
          
          UCT = Qx./Nx + c*sqrt(abs(log(Nx(Sx+1, Sy+1))./Nx)); %calculate UCT
          %restrict to four moves
            UCTi = zeros(1, 4);
                UCTi(1) = UCT(sx+1, mod1(sy+2, M)); %up
                UCTi(2) = UCT(mod1(sx+2, M), sy+1); %right 
                UCTi(3) = UCT(sx+1, mod1(sy, M)); %down
               UCTi(4) = UCT(mod1(sx, M), sy+1); %left
          [Max, R] = max(UCTi);  %Max is max, R is direction
          %move in max direction
          if R == 1
              D = up;
          elseif R == 2
              D = right; 
          elseif R == 3
              D = down; 
          else
              D = left; 
          end
          
          sx1 = sx; 
          sy1 = sy; 
          Nx(sx+1, sy+1) = Nx(sx+1, sy+1) + 1; 
          
          count = 0; 
           h = drawLevy_Hottovy(c,mu);
            h = round(h); 
            dh = 0; 
          
                while dist(Tx, Ty, sx1, sy1) > rv%random playout
                    count = count + 1; 
                    if count > Limit
                        fprintf('limit! 1\n')
                        break
                    end
                    if dh < h
                        [sx1, sy1] = movetor(sx1, sy1, D, N); 
                        dh = dh + 1; 
                    else
                    [sx1, sy1] = randmovetor(sx1, sy1, N);
                    end 
                    
              [sx1, sy1] = randmovetor(sx1, sy1, N); 
              K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %add step count
              Nx(sx1+1, sy1+1) = Nx(sx1+1, sy1+1) + 1; %add visit count to playout nodes
          end
          Kx(sx+1, sy+1) = K(sx+1, sy+1)/Nx(sx+1, sy+1); %calculate avg step count
          Qx(sx+1, sy+1) = 1/Kx(sx+1, sy+1);  %calculate win count
                    
       end
       
      
       %Play best move
            UCT = Qx./Nx + c*sqrt(abs(log(Nx(Sx+1, Sy+1))./Nx)); %calculate UCT
            %[MX, DX] = max(UCT); 
            %[MY, DY] = max(DX); 
            %BX = DY; 
            %BY = DX(DY); 
            
            %U = dist(Sx, Sy, BX, BY+1); 
            %   R = dist(Sx, Sy, BX+1, BY); 
            %   D = dist(Sx, Sy, BX, BY-1);
            %   L = dist(Sx, Sy, BX-1, BY); 
            %   O = [U, R, D, L]; 
            %   [W, A] = min(O);
            %   O = [up; right; down; left];
            %   B = O(A,:);
            %   [Sx, Sy] = movetor(Sx, Sy, B, N);
            
            
            %restrict to four moves
            UCTi = zeros(1, 4);
                UCTi(1) = UCT(Sx+1, mod1(Sy+2, M)); %up
                UCTi(2) = UCT(mod1(Sx+2, M), Sy+1); %right 
                UCTi(3) = UCT(Sx+1, mod1(Sy, M)); %down
               UCTi(4) = UCT(mod1(Sx, M), Sy+1); %left
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
    
   if length(Px) > Limit
       break
       fprintf('Limit reached')
   end
   
    end
    
    %figure
    %plot(Tx, Ty, '*r')
    %hold on
    %plot(Px, Py, 'o-b')
    %axis([0 N 0 N])
    
      MCTS(t) = length(Px);
end
Mtime = toc 



  
%% Calculate Results

MAvg = mean(MCTS); 
MVar = var(MCTS); 
MDiff = MCTS - MOPT; 
MADiff = mean(MDiff); 
MVarDiff = var(MDiff); 
MTimeAvg = Mtime/T; 

%% Print Results

fprintf('For a grid of N = %d, %d Trials, sigma = %d and %d Loops the results are as follows: \n', N, T, sigma, L)
fprintf('Average steps: %6.2f \n', MAvg)
fprintf('Variance: %7.3f \n', MVar)
fprintf('Average difference from optimal: %6.2f \n', MADiff)
fprintf('Difference variance: %7.3f \n', MVarDiff)
fprintf('Total Time: %7.4f \n', Mtime)