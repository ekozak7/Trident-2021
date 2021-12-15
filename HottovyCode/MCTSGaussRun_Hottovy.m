function [Output] = MCTSGaussRun(sigma)

%MCTS for one Target
%Searcher starts in center
%Target is known 
%Torus


% Parameters

Trials = 1;  %Number of trials
N = 40;    %Size of search grid
M = N; %Size of matrices
L = 800;  %Number of loops
c = 2;   %exploration constant
rv = 1; %vision radius
Limit = 10000; 

% Data Collection

MCTS = zeros(1, Trials);  %Number of steps per trial
MOPT = zeros(1, Trials);  %Optimal steps per trial

MAvg = mean(MCTS);   %Average steps per trial
MVar = var(MCTS);    %Variance of steps per trial


% RUN MCTS

tic
for t = 1:Trials
   t
   
   %Store Data
    Nx = zeros(M); %visit count
    Qx = zeros(M); %success rate
    K = zeros(M);  %total step count
    Kx = zeros(M); %average step count
    
    %Place Target
    Target = [gausstarget(N, sigma);gausstarget(N, sigma)];
    
    %Place Searcher
    Searcher = [1; 1];
    
    % Optimal steps
    Opt_steps = TorusDist(Searcher,Target,N);
    
    Searcher_String = Searcher;
    
    
    %Start MCTS
    while (norm(Searcher-Target)>0)
        
        %if length(Px) >= Limit
         %   fprintf('LIMIT')
          %  break
       % end
        
        % First test unvisted nodes. 
        
        % List sites
        Sites_x = [mod([Searcher(1)+1,Searcher(1)-1],N), Searcher(1),Searcher(1)]+[1,1,1,1];
        Sites_y = [Searcher(2),Searcher(2),mod([Searcher(2)+1,Searcher(2)-1],N)]+[1,1,1,1];
        
                                      
        Neighbors = [Nx(Sites_x(1),Sites_y(1)), Nx(Sites_x(2),Sites_y(2)), Nx(Sites_x(3),Sites_y(3)), Nx(Sites_x(4),Sites_y(4))];
                          
        while min(Neighbors)==0
           % Go to unvsited nodes
            [zero, Site] = min(Neighbors);
            xdir = (Site-2)*(Site-3)*(Site-4)/(-1*-2*-3)*1+(Site-1)*(Site-3)*(Site-4)/(1*-1*-2)*(-1);
            ydir = (Site-1)*(Site-2)*(Site-4)/(2*1*-1)*(1) + (Site-1)*(Site-2)*(Site-3)/(3*2*1)*(-1);

            %Set Practice Target
            Practice_Target = [gausstarget(N, sigma);gausstarget(N, sigma)];
            
            %Play out
            [Search] = Random_Playout(mod(Searcher+[xdir;ydir],N),Practice_Target,N);
            
            %LOOP AND UPDATE NODES *BackPropagate
            [New_Nx, New_Qx] = BackPropagate(Search, Nx, Qx);
            
            New_Qx(Searcher(1)+1,Searcher(2)+1) = New_Qx(Searcher(1)+1,Searcher(2)+1)-length(Search); 
            New_Nx(Searcher(1)+1,Searcher(2)+1) = New_Nx(Searcher(1)+1,Searcher(2)+1) + 1; % add visit count to root nodc
            
            %Set new sites
            Nx = New_Nx;
            Qx = New_Qx;
            
           Neighbors = [Nx(Sites_x(1),Sites_y(1)), Nx(Sites_x(2),Sites_y(2)), Nx(Sites_x(3),Sites_y(3)), Nx(Sites_x(4),Sites_y(4))];
        end
        
        %Run Simulations
        for loops = 1:L
                        
             %reset starting node 
             UCT = Qx./Nx + c*sqrt(abs(log(Nx(Searcher(1)+1,Searcher(2)+1)./Nx))); %calculate UCT
             
             Neighbors_UCT = [UCT(Sites_x(1),Sites_y(1)), UCT(Sites_x(2),Sites_y(2)), UCT(Sites_x(3),Sites_y(3)), UCT(Sites_x(4),Sites_y(4))];
            [Max, Site] = max(Neighbors_UCT);  %Max is max, R is direction
            xdir = (Site-2)*(Site-3)*(Site-4)/(-1*-2*-3)*1+(Site-1)*(Site-3)*(Site-4)/(1*-1*-2)*(-1);
            ydir = (Site-1)*(Site-2)*(Site-4)/(2*1*-1)*(1) + (Site-1)*(Site-2)*(Site-3)/(3*2*1)*(-1);
        
            % Set practice Target to Target
            Practice_Target = [gausstarget(N, sigma);gausstarget(N, sigma)];
            
            %Play out
            [Search] = Random_Playout(mod(Searcher+[xdir;ydir],N),Practice_Target,N);
            
            %LOOP AND UPDATE NODES *BackPropagate
            [New_Nx, New_Qx] = BackPropagate(Search, Nx, Qx);
            New_Qx(Searcher(1)+1,Searcher(2)+1) = New_Qx(Searcher(1)+1,Searcher(2)+1)-length(Search); 
            New_Nx(Searcher(1)+1,Searcher(2)+1) = New_Nx(Searcher(1)+1,Searcher(2)+1) + 1; % add visit count to root node
           
            %Set new sites
            Nx = New_Nx;
            Qx = New_Qx;
        end
        
        % Choose New direction
         UCT = Qx./Nx; %calculate UCT
             
        Neighbors_UCT = [UCT(Sites_x(1),Sites_y(1)), UCT(Sites_x(2),Sites_y(2)), UCT(Sites_x(3),Sites_y(3)), UCT(Sites_x(4),Sites_y(4))];
        [Max, Site] = max(Neighbors_UCT);  %Max is max, R is direction
        xdir = (Site-2)*(Site-3)*(Site-4)/(-1*-2*-3)*1+(Site-1)*(Site-3)*(Site-4)/(1*-1*-2)*(-1);
        ydir = (Site-1)*(Site-2)*(Site-4)/(2*1*-1)*(1) + (Site-1)*(Site-2)*(Site-3)/(3*2*1)*(-1);
        
        %Update Searcher
        Searcher = mod(Searcher+[xdir;ydir],N);
        
        Searcher_String = [Searcher_String Searcher];
    end
    
    Num_Steps = (length(Searcher_String)-1)
   
    MOPT(t) = Opt_steps;
    MCTS(t) = Num_Steps; 
    
end
Mtime = toc 

% Calculate Results

MAvg = mean(MCTS(isfinite(MCTS))./MOPT(isfinite(MOPT))); 
MVar = var(MCTS(isfinite(MCTS))./MOPT(isfinite(MOPT))); 
MTimeAvg = Mtime/Trials;

% Print Results

fprintf('For a grid of N = %d, %d Trials, and %d Loops the results are as follows: \n', N, Trials, L)
fprintf('Average steps (percentage): %6.2f \n', MAvg)
fprintf('Variance: %7.3f \n', MVar)
fprintf('Time per trial: %7.4f \n', MTimeAvg)
fprintf('Total time: %7.4f \n', Mtime)

% Save results

save(strcat('MCTSGaussN',num2str(N),'Loops',num2str(L),'Sigma',num2str(floor(sigma)),'Trials',num2str(Trials)))

end