% function [h] = plotf(pts,varargin) ;
% 
% Same as plot, but by taking as input matrices instead of vectors for each dimension
%
% Author: Antoine Serrurier
% Date: 25/02/2014

function [h] = plotf(PT,varargin) ;

if isempty(PT)
   h = plot([],[],varargin{:}) ; 
elseif size(PT,2) == 3
   h = plot3(PT(:,1),PT(:,2),PT(:,3),varargin{:}) ;
elseif size(PT,2) == 2
   h = plot(PT(:,1),PT(:,2),varargin{:}) ;
end

