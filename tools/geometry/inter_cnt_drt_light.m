% ALL = inter_cnt_drt_light(Xcnt, Ycnt, Xlin0, Ylin0, alpha);
% 
% Intersections between a curve (Xcnt, Ycnt) and a set of half-lines (Xlin0, Ylin0, alpha)
% 
% Inputs
%   Xcnt(1, nbPts) : X coordinates of the points of the input curve
%   Ycnt(1, nbPts) : Y coordinates of the points of the input curve
%   Xlin0(1, nbLines) : X coordinate of the origins of the lines
%   Ylin0(1, nbLines) : Y coordinate of the origins of the lines
%   alpha(1, nbLines) : Orientation of the lines in radian
% 
% Outputs
%   ALL(nbLines) : Structure containing the outputs
%       ALL(i).Intsct_X(nbPtsInter) : X coordinates of the intersection points between the curve and the line i
%       ALL(i).Intsct_Y(nbPtsInter) : Y coordinates of the intersection points between the curve and the line i
%       ALL(i).Inter(nbPtsInter)    : Distances between the intersection points and the origin point of the line i
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 24/06/2024


function TOUT = inter_cnt_drt_light(Xcnt_, Ycnt_, Xdrt0_, Ydrt0_, alpha_)

% Dimensions
nb_drt = length(alpha_);
nb_pts_cnt = length(Xcnt_);

% Constants
determ_thresh = 1.e-30;
int_thresh = 5e-10; % Threshold to defone the appartenance of an interection point to a contour

% Extremities of the segments
X1 = Xcnt_(1:nb_pts_cnt-1); Y1 = Ycnt_(1:nb_pts_cnt-1);
X2 = Xcnt_(2:nb_pts_cnt); Y2 = Ycnt_(2:nb_pts_cnt); 
% Segment contour equaltion
beta = atan2(Y2 - Y1, X2 - X1);

% Loop on the lines ----------------------------------------------------------------
for ind = 1:nb_drt
  % Equaltion of the line
  % (X-Xdrt0_)*sin(alpha_) - (Y-Ydrt0_)*cos(alpha_) = 0  ou bien
  % X*sin(alpha_) - Y*cos(alpha_) - Xdrt0_*sin(alpha_) + Ydrt0_*cos(alpha_) = 0
  
  % Possible intersections between segment and line
  Int_X = []; Int_Y = [];
  % Distances between the intersection points and the origin point of the line
  Inter_tmp = [];
  
  % Determinant
  determ = sin(alpha_(ind)) * cos(beta) - cos(alpha_(ind)) * sin(beta);
  
  ind_tmp = find(~isnan(determ));
  
  % Loop on the curve ------------------------------------------------------
  for ind_ind_tmp = 1:length(ind_tmp)
    ind2 = ind_tmp(ind_ind_tmp);
    % intersection point
    % [ind2, abs(determ(ind2)), abs(determ(ind2)) > determ_thresh]
    if abs(determ(ind2)) > determ_thresh,
        A = - Xdrt0_(ind) * sin(alpha_(ind)) + Ydrt0_(ind) * cos(alpha_(ind));
        B = - X1(ind2) * sin(beta(ind2)) + Y1(ind2) * cos(beta(ind2));
        Xinter = (- A * cos(beta(ind2)) + B * cos(alpha_(ind)))/determ(ind2);
        Yinter = (- A * sin(beta(ind2)) + B * sin(alpha_(ind)))/determ(ind2);
        % plot(Xinter, Yinter, 'cx'); drawnow
        dist_tmp = sqrt((Xinter - Xdrt0_(ind)) .^2 + (Yinter - Ydrt0_(ind)) .^2);
        % If the intersection point is between the 2 contour points...
        if ( ((X1(ind2) - Xinter) * (X2(ind2) - Xinter) <= int_thresh) & ...
                ((Y1(ind2) - Yinter) * (Y2(ind2) - Yinter) <= int_thresh)); % & ...
            Int_X = [Int_X Xinter];
            Int_Y = [Int_Y Yinter];
            Inter_tmp = [Inter_tmp dist_tmp];
        end
    else % determinant close to 0, nothing found
        Int_X = [Int_X, NaN];
        Int_Y = [Int_Y, NaN];
        Inter_tmp = [Inter_tmp, NaN];
    end
  end % for ind2 = 1:nb_pts_cnt-1 -------------------------------------------
  
  % Output
  TOUT(ind).Intsct_X = Int_X;
  TOUT(ind).Intsct_Y = Int_Y;
  TOUT(ind).Inter = Inter_tmp;
  
end % for ind = 1:nb_drt --------------------------------------------------------------



