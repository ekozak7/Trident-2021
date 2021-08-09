%MCTS for one Target
%Searcher starts in center
%Target is known 
%Torus

clear all, clc

%% Parameters

T = 1;  %Number of trials
N = 40;    %Size of search grid
M = N + 1; %Size of matrices
L = 100;  %Number of loops
c = 2;   %exploration constant
rv = 1; %vision radius
Limit = 10000; 

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
    Sx = N/2; 
    Sy = N/2;
    
    %Record Path Traveled
    Px = [Sx]; 
    Py = [Sy]; 
    
    %Measure Optimal Path
    MOPT(t) = optsteps(Tx, Ty, Sx, Sy, N); 
    
    %Start MCTS
    while dist(Sx, Sy, Tx, Ty) > rv
        if length(Px) >= Limit
            fprintf('LIMIT')
            break
        end
        
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
                while(Tx ~= sx1) || (Ty ~= sy1)%random playout
                    count = count + 1; 
                    if count > Limit
                        fprintf('limit! 1\n')
                        break
                    end
                [sx1, sy1] = randmovetor(sx1, sy1, N);
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
                while(Tx ~= sx1) || (Ty ~= sy1)  %random playout
                    count = count + 1; 
                    if count > Limit
                        fprintf('limit! 2\n')
                        break
                    end
                [sx1, sy1] = randmovetor(sx1, sy1, N); 
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
                while(Tx ~= sx1) || (Ty ~= sy1) %random playout
                    count = count + 1; 
                    if count > Limit
                        fprintf('limit! 3\n')
                        break
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
                while(Tx ~= sx1) || (Ty ~= sy1) %random playout
                    count = count + 1; 
                    if count > Limit
                        fprintf('limit! 4\n')
                        break
                    end
                [sx1, sy1] = randmovetor(sx1, sy1, N); 
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
          UCT = Qx./Nx + c*sqrt(abs(log(Nx(Sx+1, Sy+1))./Nx)); %calculate UCT
          %UCT = Qx + c*sqrt(abs(log(Nx(Sx+1, Sy+1))./Nx));
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
          while (Tx ~= sx1) || (Ty ~= sy1)
              count = count + 1; 
              if count > Limit
                  %fprintf('limit! 5\n')
                  break
              end
              [sx1, sy1] = randmovetor(sx1, sy1, N); 
              K(sx+1, sy+1) = K(sx+1, sy+1) + 1; %add step count
              Nx(sx1+1, sy1+1) = Nx(sx1+1, sy1+1) + 1; %add visit count to playout nodes
          end
          Kx(sx+1, sy+1) = K(sx+1, sy+1)/Nx(sx+1, sy+1); %calculate avg step count
          Qx(sx+1, sy+1) = 1/Kx(sx+1, sy+1);  %calculate win count
       end
       
       %Play best move
            %UCT = Qx./Nx + c*sqrt(abs(log(Nx(Sx+1, Sy+1))./Nx)); %calculate UCT
            %UCT = Qx + c*sqrt(abs(log(Nx(Sx+1, Sy+1))./Nx));
            %restrict to four moves
            %UCTi = zeros(1, 4);
            %    UCTi(1) = UCT(Sx+1, mod1(Sy+2, M)); %up
            %    UCTi(2) = UCT(mod1(Sx+2, M), Sy+1); %right 
            %    UCTi(3) = UCT(Sx+1, mod1(Sy, M)); %down
            %    UCTi(4) = UCT(mod1(Sx, M), Sy+1); %left
          %[Max, R] = max(UCTi);  %Max is max, R is direction
          %%Use max visit count
          Ni = zeros(1,4); 
            Ni(1) = Nx(Sx+1, mod1(Sy+2, M)); 
            Ni(2) = Nx(mod1(Sx+2, M), Sy+1); 
            Ni(3) = Nx(Sx+1, mod1(Sy, M)); 
            Ni(4) = Nx(mod1(Sx, M), Sy+1); 
            Ni
          [Max, R] = max(Ni); 
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
    
    MCTS(t) = length(Px); 
    end
    
end
Mtime = toc 
  
%% Calculate Results

MAvg = mean(MCTS); 
MVar = var(MCTS); 
MDiff = MCTS - MOPT; 
MADiff = mean(MDiff); 
MVarDiff = var(MDiff); 
MTimeAvg = Mtime/T;

%% Plot Example Run

%steps = length(Px); 

%figure
%plot(Tx,Ty, '*r');
%hold on
%plot(Px, Py, '-b');
%axis([0 N, 0 N])
%tit = sprintf('%d Steps', steps);
%title(tit)
%print('mctsk3', '-dpng')

%% Print Results

fprintf('For a grid of N = %d, %d Trials, and %d Loops the results are as follows: \n', N, T, L)
fprintf('Average steps: %6.2f \n', MAvg)
fprintf('Variance: %7.3f \n', MVar)
fprintf('Average difference from optimal: %6.2f \n', MADiff)
fprintf('Difference variance: %7.3f \n', MVarDiff)
fprintf('Time per trial: %7.4f \n', MTimeAvg)
fprintf('Total time: %7.4f \n', Mtime)