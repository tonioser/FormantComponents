%  [X_inter, Y_inter] = inter_drts(X_pts1, Y_pts1, X_pts2, Y_pts2, ...
%                                            X_pts3, Y_pts3, X_pts4, Y_pts4);
% 
% Returns the intersection point between 2 lines
% 
% Inputs
%   X_pts1(1): X coordinate of the first point of the first line
%   Y_pts1(1): Y coordinate of the first point of the first line
%   X_pts2(1): X coordinate of the second point of the first line
%   Y_pts2(1): Y coordinate of the second point of the first line
%   X_pts3(1): X coordinate of the first point of the second line
%   Y_pts3(1): Y coordinate of the first point of the second line
%   X_pts4(1): X coordinate of the second point of the second line
%   Y_pts4(1): Y coordinate of the second point of the second line
% 
% Outputs
%   X_inter(1): X coordinate of the intersection points
%   Y_inter(1): Y coordinate of the intersection points
%
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 24/06/2024

function [X_inter, Y_inter] = inter_drts(X_pts1, Y_pts1, X_pts2, Y_pts2, ...
    X_pts3, Y_pts3, X_pts4, Y_pts4);

determ_thresh = 1.e-30;

% Equation of the first line
% (X - X_pts1) * sin(alpha12) - (Y - Y_pts1) * cos(alpha12) = 0  or
% X * sin(alpha12) - Y * cos(alpha12) - X_pts1 * sin(alpha12) + Y_pts1 * cos(alpha12) = 0
% with 
alpha12 = atan2(Y_pts2 - Y_pts1, X_pts2 - X_pts1);

% Equation of the second line
% (X - X_pts3) * sin(alpha34) - (Y - Y_pts3) * cos(alpha34) = 0  or
% X * sin(alpha34) - Y * cos(alpha34) - X_pts3 * sin(alpha34) + Y_pts3 * cos(alpha34) = 0
% with 
alpha34 = atan2(Y_pts4 - Y_pts3, X_pts4 - X_pts3);

% Determinant
determ = sin(alpha12) * cos(alpha34) - cos(alpha12) * sin(alpha34);

if abs(determ) > determ_thresh
  A = X_pts1 * sin(alpha12) - Y_pts1 * cos(alpha12);
  B = X_pts3 * sin(alpha34) - Y_pts3 * cos(alpha34);
  X_inter = (+ A * cos(alpha34) - B * cos(alpha12)) / determ;
  Y_inter = (+ A * sin(alpha34) - B * sin(alpha12)) / determ;
else % determinant close to 0, nothing found
  X_inter = NaN;
  Y_inter = NaN;
end

