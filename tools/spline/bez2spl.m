function [Xspl, Yspl] = bez2spl(Xbez, Ybez, NbPtsSpl);
% [Xspl, Yspl] = bez2spl(Xbez, Ybez, nbPtsSpl);
% 
% Bezier function
%
% Inputs
%   Xbez(1,nbPtsCtrl) : X coordinates of the control points
%   Ybez(1,nbPtsCtrl) : Y coordinates of the control points
%   nbPtsSpl(1) : Number of points -1 for the output spline (output has nbPtsSpl+1 points)
% 
% Outputs
%   Xspl(1,nbPtsSpl+1) : X coordinates of the spline points
%   Yspl(1,nbPtsSpl+1) : Y coordinates of the spline points
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 25/06/2024


% Version issue
VERS = version;
if strcmp(VERS(1:3), '6.5')
	errordlg('bez2spl does not work with matlab 6.5!');
	error('bez2spl does not work with matlab 6.5!');
end

% At least 2 points
if length(Xbez) >= 2

  % Remove the duplicates
  ind_simple = find(diff(Xbez) ~= 0 |  diff(Ybez) ~= 0);

  % Add the last point
  ind_simple = [ind_simple, length(Xbez)];

  % Control points
  cs = cscvn_([Xbez(ind_simple); Ybez(ind_simple)]);

  % Bezier curve
  points = eval_sp(cs, NbPtsSpl); % fnplt without plot
  Xspl = points(:,1)'; Yspl = points(:,2)';

% Only one point
else
  Xspl = Xbez(1) * ones(NbPtsSpl+1,1)';
  Yspl = Ybez(1) * ones(NbPtsSpl+1,1)';
end
return
