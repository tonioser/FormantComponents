function cols = getColPlot(n)

% function cols = getColPlot(n);
% 
% Get the default colors of matlab for the desired number of objects.
% 
% Inputs
%     n(nbObj) : Number of objects
% 
% Outputs
%     cols(nbObj,3) : Matrix of colors RGB
% 
% Author: Antoine Serrurier
% Date: 10/03/2017

cols = [];

h = figure;
clf
hp = plot(rand(n,2)', rand(n,2)', '*');
for ii = 1:n
    cols = [cols; hp(ii).Color];
end  % for ii = 1:n
close(h);

return



