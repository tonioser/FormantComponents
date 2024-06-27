function cntOut = cntDilateErode(cntIn, d)

% function cntOut = cntDilateErode(cntIn, d);
% 
% Dilate or erode a 2D contour of the distance d.
% 
% The contour must be closed in order that 'dilatation' or 'erosion'
% means something. It not, the contour is artificially closed.
% 
% Inputs
%      cntIn(nbPts, iXY) : Contour
%      d(1)              : Distance to dilate (d>0) or erode (d<0) the contour
% 
% Outputs
%      cntOut(nbPts, iXY) : Dilated/Eroded contour
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 06/10/2015
%

%**************************************************
while 0
figure
FIG
plotf(cntIn);
end  % while 0
%**************************************************

% Check if the contour is closed and close it if necessary
isCntClosed = 1;
if ~isferme(cntIn(:,1)', cntIn(:,2)')
    % warning('cntDilate: the input contour is open and artificially closed for the calculations')
    cntIn = cntIn([1:end, 1], :);
    isCntClosed = 0;
end  % if ~isferme(cnt)

% Orient the contour ant-clockwise
isAntiClockwise = 1;
if polygravcentr(cntIn(:,1), cntIn(:,2)) < 0
    cntIn = cntIn(end:-1:1, :);
    isAntiClockwise = 0;
end  % if ~isferme(cnt)

% Initialisation of the output
cntOut = NaN(size(cntIn));

% For each point we dilate along its normale
for iPt = 1:size(cntIn,1)
    % Point, previous point, next point
    pt = cntIn(iPt,:);
    if iPt == 1
        pt1 = cntIn(end-1,:);
        pt2 = cntIn(2,:);
    elseif iPt == size(cntIn,1)  % if iPt == 1
        pt1 = cntIn(end-1,:);
        pt2 = cntIn(2,:);
    else  % if iPt == 1
        pt1 = cntIn(iPt-1,:);
        pt2 = cntIn(iPt+1,:);
    end  % if iPt == 1
    % Normal
    cross1 = cross([0 0 -1], norme_vecteur([pt-pt1, 0]));
    cross2 = cross([0 0 -1], norme_vecteur([pt2-pt, 0]));
    normal = (cross1(1:2) + cross2(1:2)) / 2;
    % Dilatation
    cntOut(iPt,:) = pt + d * normal;
end  % for iPt = 1:size(cnt,1)

% We re-orient it as the input
if ~isAntiClockwise
    cntOut = cntOut(end:-1:1, :);
end  % if ~isAntiClockwise

% We re-open it as the input if necessary
if ~isCntClosed
    cntOut = cntOut(1:end-1, :);
end  % if ~isCntClosed

%**************************************************
while 0
figure
FIG
plotf(cntIn);
plotf(cntOut);
end  % while 0
%**************************************************

return



