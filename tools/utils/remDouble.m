function A_unique = remDouble(A, thresh)

% function A_unique = remDouble(A {, thresh});
% 
% Remove double points in a list of points. It can be exact doubles or
% closer than a threshold to adjust.
% 
% Inputs
%      A(nbPts, iXYZ) : List of points.
% 
% Outputs
%      A_unique(nbPts, iXYZ) : List of points without doubles (same order)
% 
% Author: Antoine Serrurier
% Date: 07/07/2015

% This code removes only exact doubles
if nargin == 1
    [~, indOK] = unique(A, 'rows', 'stable');
end  % if nargin == 1

% This code removes double closer than a threshold
if nargin == 2
    nbDigits = log(1/thresh)/log(10) + 1;
    nbDigits = int8(nbDigits); % Added by AS on 06/10/2015
    [~, indOK] = unique(round(A,nbDigits), 'rows', 'stable');
end  % if nargin == 1

A_unique = A(indOK,:);

return



