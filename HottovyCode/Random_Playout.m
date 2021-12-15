function [S_final] = Random_Playout(Searcher,Target,N)
%
%
%
%
%
S = [Searcher];
S0 = Searcher;

while(norm(S0-Target)>0)

        % move
        move = randi([0 3]); % Choose random direction

        % Move spot    right (move = 1) left (move = 0)
        move_spot = (move*(2-move)*(3-move))/(1*1*2)*[1; 0] + (1-move)*(2-move)*(3-move)/(1*2*3)*[-1; 0] + move*(1-move)*(3-move)/(2*-1*1)*[0; 1] + move*(1-move)*(2-move)/(3*-2*-1)*[0; -1];

        S0 = mod(S0+move_spot,N);
        S = [S S0];

end

S_final = S;