function [cntNomos, predNomos, indMean] = nomos2D(moy, pred, coef, nbSteps, linORsin, boolMean)

% function [cntNomos, scoreNomos, indMean] = nomos2D(meanCnts, scores, basisVectors {, nbSteps, linORsin, boolMean});
% 
% Calculation of 2D contour nomograms in linear models.
% 
% The nomograms are calculated between the minimal and maximal values of
% the predictors.
% 
% Inputs
%      meanCnts(nbPts, iXY)             : Mean contour
%      scores(nbObs, nbPred)            : Scores
%      basisVectors(nbPred, nbPts, iXY) : Basis vectors
%                                         if nbPred=1, then the size can be basisVectors(nbPts, iXY)
%      nbSteps(int)                     : Number of steps of the nomograms [default = 10] (skipped by using NaN, or Inf or []) 
%      linORsin(string)                 : Nomograms on linar stepas ('lin') or sinusoidal ('sin') [default = 'lin'] (skipped by using NaN, or Inf or []) 
%      boolMean(1)                      : True (=1) if we want also to include the case where score=0 [default = 0] (skipped by using NaN, or Inf or [])
% 
% Outputs
%      cntNomos(nbPred, nbSteps, nbPts, iXY) : Nomogram contours 
%      scoreNomos(nbPred, nbSteps)           : Score values used for the nomograms
%      indMean(nbSteps)                      : Index in scoreNomos corresponding to the nomogram calculated with score=0 (calculated when boolMean==1) 
% 
% Author: Antoine Serrurier
% Date (adaptation): 25/06/2024

nbSteps_default = 10;
linORsin_default = 'lin';
boolMean_default = 0;
if nargin < 4
    nbSteps = nbSteps_default;
end  % if nargin < 4
if isnan(nbSteps) | isinf(nbSteps) | isempty(nbSteps)
    nbSteps = nbSteps_default;
end  % if isnan(nbSteps) | isinf(nbSteps) | isempty(nbSteps)
if nargin < 5
    linORsin = linORsin_default;
end  % if nargin < 4
if isnan(linORsin) | isinf(linORsin) | isempty(linORsin)
    linORsin = linORsin_default;
end  % if isnan(nbSteps) | isinf(nbSteps) | isempty(nbSteps)
if nargin < 6
    boolMean = boolMean_default;
end  % if nargin < 4
if isnan(boolMean) | isinf(boolMean) | isempty(boolMean)
    boolMean = boolMean_default;
end  % if isnan(nbSteps) | isinf(nbSteps) | isempty(nbSteps)

% Case where there is only 1 dimension
[dim1, dim2] = size(moy);
if (dim1 == 1) || (dim2 == 1)
    dim = 1; % 1D nomograms
else  % if (dim1 == 1) || (dim2 == 1)
    % dim = 2; % 2D nomograms
    dim = dim2; % 2D nomograms % Written on 02/08/2016
end  % if (dim1 == 1) || (dim2 == 1)

% Case where there is only one score (= 1 predictor)
if dim == 1
    [dim1, dim2] = size(coef);
    if (dim1 == 1) || (dim2 == 1)
        % It must be a line vector
        if dim2 == 1
            coef = coef';
        end  % if dim2 == 1
    end  % if (dim1 == 1) || (dim2 == 1)
end  % if dim == 1
% if dim == 2
if dim >= 2 % Written on 02/08/2016
    if length(size(coef)) == 2
        coef = reshape(coef, [1 size(coef)]);
    end  % if length(size(coef)) == 2
end  % if dim == 2

% Dimensions
nbObs = size(pred, 1);
nbPred = size(pred, 2);
nbPts = size(coef, 2);

% Init output
cntNomos = NaN([nbPred, nbSteps, nbPts, dim]);
predNomos = NaN(nbPred, nbSteps);
if boolMean
    cntNomos = NaN([nbPred, nbSteps+1, nbPts, dim]);
    predNomos = NaN(nbPred, nbSteps+1);
end  % if boolMean
indMean = NaN(nbPred);

% Loop
for iPred = 1:nbPred
    % Predictor values
    predMin = min(pred(:,iPred));
    predMax = max(pred(:,iPred));
    switch linORsin
        case 'lin'  % switch linORsin
            predVal = linspace(predMin,predMax,nbSteps);
        case 'sin'  % switch linORsin
            predVal = sin(linspace(-pi/2,pi/2,nbSteps)) * (predMax - predMin) / 2 + mean([predMin, predMax]);
    end  % switch linORsin
    if boolMean
        [predVal, indSort] = sort([predVal, 0]);
        indMean(iPred) = find(indSort == (nbSteps+1));
    end  % if boolMean
    
    % Prediction
    for iStep = 1:length(predVal)
        cntNomos(iPred, iStep, :, :) = moy + predVal(iStep) * squeeze(coef(iPred,:,:));
    end  % for iStep = 1:nbSteps
    % Predictor values
    predNomos(iPred, :) = predVal;
end  % for iPred = 1:5


return



