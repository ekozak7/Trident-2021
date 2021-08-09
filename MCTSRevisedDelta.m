%Revised MCTS Code
%Using different final selection policies
%Delta Distribution
%7/19/21

clear all, clc

%% Parameters 

T = 1; 
N = 40; 
M = N+1; 
L = 100; 
c = 2; 
rv = 1; 
TLimit = 10000; 

up = [0 1]; 
right = [1 0]; 
down = [0 -1];
left = [-1 0];

%% Data 

MCTS = zeros(1,T); %Steps per trial
MOPT = zeros(1,T); 

%% Run Trials
 
tic
for t = 1:T
 t
 
 %Store Data
 Nx = zeros(M); %Visit count
 Qx = zeros(M); %Reward
 K = ones(M); %Step count
 
 %Target
 Tx = N/2;
 Ty = N/2; 
 %Searcher
 Sx = 1; 
 Sy = 1;  
 %Path Traveled
 Px = [Sx]; 
 Py = [Sy]; 
 %Optimal Path
 MOPT(t) = optsteps(Tx, Ty, Sx, Sy, N); 
 
 %Start MCTS
 while dist(Sx, Sy, Tx, Ty) > rv
      Nx = zeros(M); %Visit count
      Qx = zeros(M); %Reward
      K = ones(M); %Step count
      
     %Unvisited Nodes
     if Nx(Sx+1, mod1(Sy+2, M)) == 0  %up
         Nx(Sx+1, Sy+1) = Nx(Sx+1, Sy+1) + 1;
         [sx, sy] = movetor(Sx, Sy, up, N);
         Nx(sx+1, sy+1) = Nx(sx+1, sy+1) + 1; 
         sx1 = sx; 
         sy1 = sy; 
         tx1 = N/2; 
         ty1 = N/2;
         count = 0; 
         while dist(tx1, ty1, sx1, sy1) > rv
             count = count + 1; 
             if count>TLimit
                 break
             end
             [sx1, sy1] = randmovetor(sx1, sy1, N); 
             %K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %step count of selected node
         end
         K(sx+1, sy+1) = K(sx+1, sy+1) + 1/count;
     end
       
      if Nx(mod1(Sx+2,M), Sy+1) == 0  %right
          Nx(Sx+1, Sy+1) = Nx(Sx+1, Sy+1) + 1;
         [sx, sy] = movetor(Sx, Sy, right, N);
         Nx(sx+1, sy+1) = Nx(sx+1, sy+1) + 1; 
         sx1 = sx; 
         sy1 = sy; 
         tx1 = N/2; 
         ty1 = N/2;
         count = 0; 
         while dist(tx1, ty1, sx1, sy1) > rv
             count = count + 1; 
             if count>TLimit
                 break
             end
             [sx1, sy1] = randmovetor(sx1, sy1, N); 
             %K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %step count of selected node
         end 
         K(sx+1, sy+1) = K(sx+1, sy+1) + 1/count;
      end
     
      if Nx(Sx+1, mod1(Sy, M)) == 0  %down
          Nx(Sx+1, Sy+1) = Nx(Sx+1, Sy+1) + 1;
         [sx, sy] = movetor(Sx, Sy, down, N);
         Nx(sx+1, sy+1) = Nx(sx+1, sy+1) + 1; 
         sx1 = sx; 
         sy1 = sy; 
         tx1 = N/2; 
         ty1 = N/2;
         count = 0; 
         while dist(tx1, ty1, sx1, sy1) > rv
             count = count + 1; 
             if count>TLimit
                 break
             end
             [sx1, sy1] = randmovetor(sx1, sy1, N); 
             %K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %step count of selected node
         end  
         K(sx+1, sy+1) = K(sx+1, sy+1) + 1/count;
      end
         
      if Nx(mod1(Sx,M), Sy+1) == 0  %left
          Nx(Sx+1, Sy+1) = Nx(Sx+1, Sy+1) + 1;
         [sx, sy] = movetor(Sx, Sy, left, N);
         Nx(sx+1, sy+1) = Nx(sx+1, sy+1) + 1; 
         sx1 = sx; 
         sy1 = sy; 
         tx1 = N/2; 
         ty1 = N/2;
         count = 0; 
         while dist(tx1, ty1, sx1, sy1) > rv
             count = count + 1; 
             if count>TLimit
                 break
             end
             [sx1, sy1] = randmovetor(sx1, sy1, N); 
             %K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %step count of selected node
         end  
         K(sx+1, sy+1) = K(sx+1, sy+1) + 1/count;
      end
     
      %Run simulations
      for loops = 1:L
          sx = Sx; 
          sy = Sy; 
          Nx(Sx+1, Sy+1) = Nx(Sx+1, Sy+1) + 1;
          %Qx = 1./K; %Total reward
          Qx = K; 
          UCT = Qx./Nx + c*sqrt(abs(log(Nx(Sx+1, Sy+1))./Nx));
          UCTi = zeros(1, 4); 
              UCTi(1) = UCT(sx+1, mod1(sy+2,M)); %up
              UCTi(2) = UCT(mod1(sx+2, M), sy+1); %right 
              UCTi(3) = UCT(sx+1, mod1(sy, M)); %down
              UCTi(4) = UCT(mod1(sx, M), sy+1); %left
          [Max, R] = max(UCTi);
          if R == 1
              [sx,sy] = movetor(sx,sy, up, N); 
          elseif R == 2
              [sx,sy] = movetor(sx,sy, right, N);
          elseif R == 3
              [sx,sy] = movetor(sx,sy, down, N);
          else
              [sx,sy] = movetor(sx,sy, left, N);
          end
          Nx(sx+1, sy+1) = Nx(sx+1, sy+1) + 1;
          sx1 = sx; 
          sy1 = sy;
          %random playout
          tx1 = N/2; 
         ty1 = N/2;
         count = 0; 
         while dist(tx1, ty1, sx1, sy1) > rv
             count = count + 1; 
             if count>TLimit
                 break
             end
             [sx1, sy1] = randmovetor(sx1, sy1, N); 
             %K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %step count of selected node
         end 
         K(sx+1, sy+1) = K(sx+1, sy+1) + 1/count;
         
      end %end loops
      
      
     %Play best move
        %Use Max Reward
        %Qx = 1./K;  %Total reward
        QN = K./Nx; %Average reward
        Qi(1) = QN(Sx+1, mod1(Sy+2, M)); 
        Qi(2) = QN(mod1(Sx+2, M), Sy+1); 
        Qi(3) = QN(Sx+1, mod1(Sy, M)); 
        Qi(4) = QN(mod1(Sx, M), Sy+1); 
        [Max, R] = max(Qi);  
        
        %Use max visit count
        %Ni(1) = Nx(Sx+1, mod1(Sy+2, M)); 
        %Ni(2) = Nx(mod1(Sx+2, M), Sy+1); 
        %Ni(3) = Nx(Sx+1, mod1(Sy, M)); 
        %Ni(4) = Nx(mod1(Sx, M), Sy+1); 
        %[Max, R] = max(Ni); 
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

 end

 MCTS(t) = length(Px); 
end

MTime = toc; 

%% Calculate Results

MDiff = MCTS - MOPT; 
MAvg = mean(MCTS); 
MADiff = mean(MDiff); 
MVar = var(MDiff); 
MTimeAvg = MTime/T; 

%% Print Results
fprintf('The average steps is: %12.5f\n', MAvg)
fprintf('The average difference over optimal is: %12.5f\n', MADiff)
fprintf('The variance is: %12.5f\n', MVar)
fprintf('The average time per trial is: %12.5f\n', MTimeAvg)

