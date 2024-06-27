function D = distdiff(V)

% function D = distdiff(V);
% 
% Distance of each segment in a vector V
% 
% Inputs
%      V(nbPts, iXYZ) : Vector of points.
% 
% Outputs
%      D(nbPts-1) : Segment distances 
% 
% Author: Antoine Serrurier
% Date: 22/07/2015

D = sqrt(sum((V(2:end,:) - V(1:end-1,:)).^2,2));

return



