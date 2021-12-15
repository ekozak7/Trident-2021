function [New_Nx,New_Qx] = BackPropagate(Search_String, Old_Nx, Old_Qx)
%
%
%   Function to back propagate results and give a new reward function for
%   each site (Qx) and new site visit function for each site (Nx)
%
[~,num_steps] = size(Search_String);

% Loop through search string. Note that the first entry is the searcher 
% position and the last entry of search string is the target site. 

for i = 1:(num_steps)
    % increment the Simulation Count
   Old_Nx(Search_String(1,i)+1,Search_String(2,i)+1)  = Old_Nx(Search_String(1,i)+1,Search_String(2,i)+1)+1;
   % Set the reward funciton
   
   Old_Qx(Search_String(1,i)+1,Search_String(2,i)+1) = Old_Qx(Search_String(1,i)+1,Search_String(2,i)+1) - (num_steps-(i));
    
end

New_Nx = Old_Nx;
New_Qx = Old_Qx;

end
    
    