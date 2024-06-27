function Pts = intersection_line_contour(Line, C)
% 
% function Pts = intersection_line_contour(Line, C);
%
% Intersection between a line and a curve 2D.
% 
% Inputs
%    Line(1:2,iXY) : 2D line
%           Line(1,iXY) : point on the line
%           Line(2,iXY) : direction vector
%    C(nb_pts,iXY)   : 2D curve
% 
% Outputs
%    Pts(nbPts,iXY) : intersection points
% 
% Author : Antoine Serrurier
% Date: 21/06/2024

% Data format
Xcnt_ = C(:,1)';
Ycnt_ = C(:,2)';
Xdrt0_ = Line(1,1);
Ydrt0_ = Line(1,2);
alpha_ = atan2(Line(2,2),Line(2,1));
sens_ = ones(size(Xdrt0_));


% Initialisations
nb_drt = length(alpha_);
nb_pts_cnt = length(Xcnt_);

% Constants
determ_thresh = 1.e-30;
int_thresh = 5e-10; % Threshold to determine if a point belongs to a contour

% Extremities
X1 = Xcnt_(1:nb_pts_cnt-1); Y1 = Ycnt_(1:nb_pts_cnt-1);
X2 = Xcnt_(2:nb_pts_cnt); Y2 = Ycnt_(2:nb_pts_cnt); 
% Equation of the region
beta = atan2(Y2 - Y1, X2 - X1);

% Loops on the lines ----------------------------------------------------------------
% Only one line here...
for ind = 1:nb_drt
  
  % possible intersections between the contour and the grid
  Int_X = []; Int_Y = [];
  % Distance between the reference point on the grid and the contour point
  Inter_tmp = [];
  
  % Determinant
  determ = sin(alpha_(ind)) * cos(beta) - cos(alpha_(ind)) * sin(beta);
  
  ind_tmp = find(~isnan(determ));
  
  % loop on the contour ------------------------------------------------------
  for ind_ind_tmp = 1:length(ind_tmp)
    ind2 = ind_tmp(ind_ind_tmp);
    % Intersection point is:
    % [ind2, abs(determ(ind2)), abs(determ(ind2)) > determ_thresh]
    if abs(determ(ind2)) > determ_thresh,
        A = - Xdrt0_(ind) * sin(alpha_(ind)) + Ydrt0_(ind) * cos(alpha_(ind));
        B = - X1(ind2) * sin(beta(ind2)) + Y1(ind2) * cos(beta(ind2));
        Xinter = (- A * cos(beta(ind2)) + B * cos(alpha_(ind)))/determ(ind2);
        Yinter = (- A * sin(beta(ind2)) + B * sin(alpha_(ind)))/determ(ind2);
        % plot(Xinter, Yinter, 'cx'); drawnow
        dist_tmp = sqrt((Xinter - Xdrt0_(ind)) .^2 + (Yinter - Ydrt0_(ind)) .^2);
        % if the intersection point is between the two contour points ...
        if ( ((X1(ind2) - Xinter) * (X2(ind2) - Xinter) <= int_thresh) & ...
                ((Y1(ind2) - Yinter) * (Y2(ind2) - Yinter) <= int_thresh)); % & ...
            Int_X = [Int_X Xinter];
            Int_Y = [Int_Y Yinter];
            Inter_tmp = [Inter_tmp dist_tmp];
        end
    else % determinant close to 0, no point found
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

Pts = [TOUT.Intsct_X', TOUT.Intsct_Y'];

end

