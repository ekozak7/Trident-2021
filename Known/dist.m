function [d] = dist(tx, ty, sx, sy)
%computes distance from searcher to target
d = sqrt((tx - sx)^2 + (ty - sy)^2);

end

