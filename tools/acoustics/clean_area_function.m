function [area_clean, length_clean] = clean_area_function(area, length, max_length_tube)

%
% function [area_clean, length_clean] = clean_area_function(area, length {, max_length_tube}); 
% 
% Split into the area function the tubes longer than max_length_tube
% 
% Input
%	area(nb_tubes)
%	length(nb_tubes)
%	max_length_tube (default: [0.5])
%
% Outputs
%	area_clean(nb_tubes)
%	length_clean(nb_tubes)
%
% Author: Antoine Serrurier
% Date: 28/04/2008

if nargin < 3
	max_length_tube = 0.5;
end  % if nargin < 3

% Vectors in lines
[dim1, dim2] = size(area);
if dim1 > dim2
	area = area';
	length = length';
end  % if dim1 > dim2

% Find tubes longer than max_length_tube
ind_tubes_long = find(length>max_length_tube);

while ~isempty(ind_tubes_long)
	i_tube = ind_tubes_long(1);
	
	% Number of new tubes
	nb_new_tubes = ceil(length(i_tube) / max_length_tube);
	
	% New area function
	area = [area(1:i_tube-1), repmat(area(i_tube),1,nb_new_tubes), area(i_tube+1:end)];
	length = [length(1:i_tube-1), repmat(length(i_tube)/nb_new_tubes,1,nb_new_tubes), length(i_tube+1:end)];
	
	% Is there any other tube to split?
	ind_tubes_long = find(length>max_length_tube);
	
end  % while ~isempty(ind_tubes_long)

% Output
area_clean = area;
length_clean = length;

if dim1 > dim2
	area_clean = area_clean';
	length_clean = length_clean';
end  % if dim1 > dim2

return
