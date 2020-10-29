%MCTS for one Target
%Searcher starts in center
%Target is known 
%Walled scenario (not a torus)

%% Parameters

T = 100;  %Number of trials
N = 40;    %Size of search grid
L = 500;  %Number of loops
c = 2;   %exploration constant

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
    Nx = zeros(N+1); %visit count
    Qx = zeros(N+1); %success rate
    K = zeros(N+1);  %total step count
    Kx = zeros(N+1); %average step count
    
    %Place Target
    [Tx, Ty] = randtarget2(N); 
    
    %Place Searcher
    Sx = N/2; 
    Sy = N/2;
    
    %Record Path Traveled
    Px = [Sx]; 
    Py = [Sy]; 
    
    %Measure Optimal Path
    MOPT(t) = optsteps(Tx, Ty, Sx, Sy, N); 
    
    %Start MCTS
    while (Sx~= Tx) || (Sy ~= Ty)
       
        %%Test unvisited nodes
        %First node: UP
        if Sy < N
            if Nx(Sx, Sy + 1) == 0  
            [sx, sy] = move1(Sx, Sy, 1, N);  %move up
            Nx(Sx, Sy) = Nx(Sx, Sy) + 1;  %add visit count to root node
            Nx(Sx, Sy+1) = Nx(Sx, Sy+1) + 1; %add visit count to up node
            sx1 = sx; %create disposable variables
            sy1 = sy; 
                while(Tx ~= sx1) || (Ty ~= sy1) %random playout
                [sx1, sy1] = randmove1(sx1, sy1, N); 
                K(sx, sy) = K(sx, sy) + 1; %add step count to selected node
                Nx(sx1, sy1) = Nx(sx1, sy1) + 1; %add visit count to playout nodes
                end
            Kx(sx, sy) = K(sx, sy)/Nx(sx, sy);  %calculate avg step count
            Qx(sx, sy) = 1/Kx(sx, sy);     %update win count
            end
        end
        %Second node: Right
        if Sx < N 
            if Nx(Sx + 1, Sy) == 0  
            [sx, sy] = move1(Sx, Sy, 2, N);  %move up
            Nx(Sx, Sy) = Nx(Sx, Sy) + 1;  %add visit count to root node
            Nx(Sx + 1, Sy) = Nx(Sx + 1, Sy) + 1; %add visit count to up node
            sx1 = sx; %create disposable variables
            sy1 = sy; 
                while(Tx ~= sx1) || (Ty ~= sy1) %random playout
                [sx1, sy1] = randmove1(sx1, sy1, N); 
                K(sx, sy) = K(sx, sy) + 1; %add step count to selected node
                Nx(sx1, sy1) = Nx(sx1, sy1) + 1; %add visit count to playout nodes
                end
            Kx(sx, sy) = K(sx, sy)/Nx(sx, sy);  %calculate avg step count
            Qx(sx, sy) = 1/Kx(sx, sy);     %update win count
            end
        end
        %Third node: Down
        if Sy >1 
            if Nx(Sx, Sy - 1) == 0  
            [sx, sy] = move1(Sx, Sy, 3, N);  %move down
            Nx(Sx, Sy) = Nx(Sx, Sy) + 1;  %add visit count to root node
            Nx(Sx, Sy - 1) = Nx(Sx, Sy - 1) + 1; %add visit count to up node
            sx1 = sx; %create disposable variables
            sy1 = sy; 
                while(Tx ~= sx1) || (Ty ~= sy1) %random playout
                [sx1, sy1] = randmove1(sx1, sy1, N); 
                K(sx, sy) = K(sx, sy) + 1; %add step count to selected node
                Nx(sx1, sy1) = Nx(sx1, sy1) + 1; %add visit count to playout nodes
                end
            Kx(sx, sy) = K(sx, sy)/Nx(sx, sy);  %calculate avg step count
            Qx(sx, sy) = 1/Kx(sx, sy);     %update win count
            end
        end
        %Fourth node: Left
        if Sx > 1
            if Nx(Sx-1, Sy) == 0  
            [sx, sy] = move1(Sx, Sy, 4, N);  %move left
            Nx(Sx, Sy) = Nx(Sx, Sy) + 1;  %add visit count to root node
            Nx(Sx-1 , Sy) = Nx(Sx-1, Sy) + 1; %add visit count to up node
            sx1 = sx; %create disposable variables
            sy1 = sy; 
                while(Tx ~= sx1) || (Ty ~= sy1) %random playout
                [sx1, sy1] = randmove1(sx1, sy1, N); 
                K(sx, sy) = K(sx, sy) + 1; %add step count to selected node
                Nx(sx1, sy1) = Nx(sx1, sy1) + 1; %add visit count to playout nodes
                end
            Kx(sx, sy) = K(sx, sy)/Nx(sx, sy);  %calculate avg step count
            Qx(sx, sy) = 1/Kx(sx, sy);     %update win count
            end
        end
        
       %Run simulations
       for loops = 1:L
          sx = Sx;  %reset starting node 
          sy = Sy;
          Nx(sx, sy) = Nx(sx, sy) + 1;  %add visit count to root node
          UCT = Qx./Nx + c*sqrt(abs(log(Nx(Sx, Sy))./Nx)); %calculate UCT
          %restrict to four movements
          UCTi = zeros(1, 4); 
          if (sy + 1) > N  %up not possible
              UCTi(1) = 0;
          else
              UCTi(1) = UCT(sx, sy+1); %up is possible
          end
          if (sx + 1) > N %right not possible
              UCTi(2) = 0; 
          else
              UCTi(2) = UCT(sx+1, sy); %right is possible
          end
          if (sy - 1) < 1 %down not possible
              UCTi(3) = 0; 
          else
              UCTi(3) = UCT(sx, sy-1); %down is possible
          end
          if (sx-1) < 1 %left not possible
              UCTi(4) = 0; 
          else 
              UCTi(4) = UCT(sx-1, sy); %left is possible
          end
          [M, R] = max(UCTi);  %M is max, R is direction
          %move in max direction
          [sx, sy] = move1(sx, sy, R, N);
          %disposable variables
          sx1 = sx; 
          sy1 = sy; 
          Nx(sx, sy) = Nx(sx, sy) + 1;  %add visit count for selected node
          %random playout
          while (Tx ~= sx1) || (Ty ~= sy1)
              [sx1, sy1] = randmove1(sx1, sy1, N); 
              K(sx, sy) = K(sx, sy) + 1; %add step count
              Nx(sx1, sy1) = Nx(sx1, sy1) + 1; %add visit count to playout nodes
          end
          Kx(sx, sy) = K(sx, sy)/Nx(sx, sy); %calculate avg step count
          Qx(sx, sy) = 1/Kx(sx, sy);  %calculate win count
       end
       
       %Play best move
            UCT = Qx./Nx + c*sqrt(abs(log(Nx(Sx, Sy))./Nx)); %calculate UCT
            %restrict to four moves
            UCTi = zeros(1, 4); 
          if (Sy + 1) > N  %up not possible
              UCTi(1) = 0;
          else
              UCTi(1) = UCT(Sx, Sy+1); %up is possible
          end
          if (Sx + 1) > N %right not possible
              UCTi(2) = 0; 
          else
              UCTi(2) = UCT(Sx+1, Sy); %right is possible
          end
          if (Sy - 1) < 1 %down not possible
              UCTi(3) = 0; 
          else
              UCTi(3) = UCT(Sx, Sy-1); %down is possible
          end
          if (Sx-1) < 1 %left not possible
              UCTi(4) = 0; 
          else 
              UCTi(4) = UCT(Sx-1, Sy); %left is possible
          end
          %Maximize UCT
          [M, R] = max(UCTi);  %M is max, R is direction
          %move in max direction
          [Sx, Sy] = move1(Sx, Sy, R, N);
          %Add to path
          Px = [Px Sx];
          Py = [Py Sy]; 
    end
    
    MCTS(t) = length(Px);
    
    %if t <= 5
       % mx = [0:20]; 
       % my = mx; 
       %[X, Y] = meshgrid(mx, my); 
       %Z = UCT;
       %steps = MCTS(t);
       %tit = sprintf('Target at x = %d, y = %d with %d steps', Tx, Ty, steps)
       %figure
       %surf(X,Y,Z)
       %title(tit)
    %end
           
end
Mtime = toc; 
  
%% Calculate Results

MAvg = mean(MCTS); 
MVar = var(MCTS); 
MDiff = MCTS - MOPT; 
MADiff = mean(MDiff); 
MVarDiff = var(MDiff); 
MTimeAvg = Mtime/T; 

%% Print Results

fprintf('For a grid of N = %d, %d Trials, and %d Loops the results are as follows: \n', N, T, L)
fprintf('Average steps: %6.2f \n', MAvg)
fprintf('Variance: %7.3f \n', MVar)
fprintf('Average difference from optimal: %6.2f \n', MADiff)
fprintf('Difference variance: %7.3f \n', MVarDiff)
fprintf('Time per trial: %7.4f \n', MTimeAvg)



