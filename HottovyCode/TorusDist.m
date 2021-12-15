function [Norm] = TorusDist(X,Y,N); 

Norm = min([norm(X-Y,1),N-2*abs(X(1)-Y(1))+norm(X-Y,1),N-2*abs(X(2)-Y(2))+norm(X-Y,1),2*N-norm(X-Y,1)]);

end