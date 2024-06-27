function [S, Xgc, Ygc] = polygravcentr(x, y, mode)
% [S, Xgc, Ygc] = polygravcentr(x, y{, mode});
% 
% Returns the signed or unsigned surface and gravity centres of curves.
% Positive if the curve is clockwise, negative otherwise.
%
% Note:
%   - Curves should not be closed
%   - Each curve should not intersect itself
% 
% Triangulation with commmon origin of the triangles at (0,0)
% Equivalent to the integral Riemann formula for surfaces
% 
% Inputs:
%   x(nbPts,nbCurves) : X coordinates of the input curves
%   y(nbPts,nbCurves) : Y coordinates of the input curves
%   mode(str)         : Signed ('sign') or unsigned ('abs') surface [default = 'abs']
% 
% Outputs:
%   S(1,nbCurves)   : Surfaces
%   Xgc(1,nbCurves) : X coordinates of the gravity centres of all curves
%   Ygc(1,nbCurves) : Y coordinates of the gravity centres of all curves
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 24/06/2024

if ~(size(x) == size(y)), error('X and Y must have the same size.'); end
if isempty(x) | isempty(y)
  S=0; Xgc = NaN; Ygc = NaN;
  return
end

% Triangle surfaces [Orig (0,0); Pt_k, Pt_kp1]
Sk = ( (x .* y([2:size(x,1) 1],:)) - ...
    (x([2:size(y,1) 1],:) .* y) ) /2;

% Curve surfaces
S = sum(Sk);

% Triangle gravity centres [Orig (0,0); Pt_k, Pt_kp1]
Xgk = (x + x([2:size(x,1) 1],:))/3;
Ygk = (y + y([2:size(y,1) 1],:))/3;

% Curve gravity centres
if S == 0
	Xgc = sum(Xgk .* Sk) .* Inf;
	Ygc = sum(Ygk .* Sk) .* Inf;
else  % if S == 0
	Xgc = sum(Xgk .* Sk) ./ S;
	Ygc = sum(Ygk .* Sk) ./ S;
end  % if S == 0

if nargin > 2
	if strcmp(mode, 'abs')
	  Xgc(find(isnan(Xgc))) = 0;
	  Ygc(find(isnan(Ygc))) = 0;
	  S(find(isnan(S))) = 0;
	  S = abs(S);
	end
end
return




