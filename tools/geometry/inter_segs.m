%  [X_inter, Y_inter] = inter_segs(X_pts1, Y_pts1, X_pts2, Y_pts2, ...
%                                            X_pts3, Y_pts3, X_pts4, Y_pts4);
% 
% Returns the intersection point between 2 segments
% 
% Inputs
%   X_pts1(1): X coordinate of the first point of the first segment
%   Y_pts1(1): Y coordinate of the first point of the first segment
%   X_pts2(1): X coordinate of the second point of the first segment
%   Y_pts2(1): Y coordinate of the second point of the first segment
%   X_pts3(1): X coordinate of the first point of the second segment
%   Y_pts3(1): Y coordinate of the first point of the second segment
%   X_pts4(1): X coordinate of the second point of the second segment
%   Y_pts4(1): Y coordinate of the second point of the second segment
% 
% Outputs
%   X_inter(1): X coordinate of the intersection points
%   Y_inter(1): Y coordinate of the intersection points
%
% Author: Antoine Serrurier
% Date (adaptation): 24/06/2024

function [X_inter, Y_inter] = inter_segs(X_pts1, Y_pts1, X_pts2, Y_pts2, ...
    X_pts3, Y_pts3, X_pts4, Y_pts4);


% Line intersections
[X_inter, Y_inter] = inter_drts(X_pts1, Y_pts1, X_pts2, Y_pts2, ...
    X_pts3, Y_pts3, X_pts4, Y_pts4);

if ~isnan(X_inter)
	
	% Verify that the intersection points belong to the segment
	% Segment 1
	vect1 = [X_inter - X_pts1; Y_inter - Y_pts1];
	vect2 = [X_inter - X_pts2; Y_inter - Y_pts2];
	% Segment 2
	vect3 = [X_inter - X_pts3; Y_inter - Y_pts3];
	vect4 = [X_inter - X_pts4; Y_inter - Y_pts4];
	
	if (dot(vect1,vect2) >= 0) || (dot(vect3,vect4) >= 0)
		X_inter = NaN; Y_inter = NaN;
	end  % if (dot(vect1,vect2) >= 0) || (dot(vect3,vect4) >= 0)
	
end  % if ~isnan(X_inter)




