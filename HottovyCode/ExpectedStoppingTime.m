%% Code to run random walk on 2D lattice to see confidence interval for expected stopping time

N = 20; % Lattice size

Total_Trial = 800;

Time = zeros(1,Total_Trial);
tic
for trial = 1:Total_Trial
   
    S0 = [5;6]; % initial conditions
    S = [S0];
    % Choose target
    T = [2;2];


    while(norm(S0-T)>0)

        % move
        move = randi([0 3]); % Choose random direction

        % Move spot    right (move = 1) left (move = 0)
        move_spot = (move*(2-move)*(3-move))/2*[1; 0] + (1-move)*(2-move)*(3-move)/6*[-1; 0] + move*(1-move)*(3-move)/(-2)*[0; 1] + move*(1-move)*(2-move)/(-6)*[0; -1];

        S0 = mod(S0+move_spot,N);
        S = [S S0];

    end

    Time(trial) = length(S);
   
end
toc
%%

SEM = std(Time)/sqrt(length(Time));
ts = tinv([0.025  0.975],length(Time)-1);
CI = mean(Time) + ts*SEM
   
Time = zeros(1,Total_Trial);
tic
for trial = 1:Total_Trial
   
    S0 = [5;4]; % initial conditions
    S = [S0];
    % Choose target
    T = [2; 2];


    while(norm(S0-T)>0)

        % move
        move = randi([0 3]); % Choose random direction

        % Move spot    right (move = 1) left (move = 0)
        move_spot = (move*(2-move)*(3-move))/2*[1; 0] + (1-move)*(2-move)*(3-move)/6*[-1; 0] + move*(1-move)*(3-move)/(-2)*[0; 1] + move*(1-move)*(2-move)/(-6)*[0; -1];

        S0 = mod(S0+move_spot,N);
        S = [S S0];

    end

    Time(trial) = length(S);
   
end
toc
%%

SEM = std(Time)/sqrt(length(Time));
ts = tinv([0.025  0.975],length(Time)-1);
CI = mean(Time) + ts*SEM
   