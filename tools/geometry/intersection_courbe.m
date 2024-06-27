function [Pts,Ou] = intersection_courbe(C1,C2) ;
%
% Returns the intersection points between two 2D curves.
% 
% If only one curve is provided, it returns the points where the curve crossed itself
% 
% Inputs
%   C1(nbPts1, iXY): Coordinates of the points of the first curve
%   C2(nbPts2, iXY): Coordinates of the points of the second curve (optional)
% 
% Outputs
%   Pts(nbPtInter, iXY): Coordinates of the intersection points
%   Ou : Localisation of the intersection points in the curves
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 24/06/2024

% Input argulents
if nargin == 1 ;
    % cas d'une seule courbe
    C2 = C1 ;
end

% Intersections

% a) Data formatting
dC1 = diff(C1,1,1) ; % M1M2 vectors
dC2 = diff(C2,1,1) ; % PiPi+1 vectors

% b) Intersection matrix 
% ---> Segment M1M2
B = dC1(:,1) * dC2(:,2)' - dC1(:,2) * dC2(:,1)' ;
A = C1(1:end-1,1) * dC2(:,2)' - C1(1:end-1,2) * dC2(:,1)' - ...
    ones(size(dC1,1),1) * (C2(1:end-1,1) .* C2(2:end,2) - C2(1:end-1,2) .* C2(2:end,1))' ;
A(find(abs(A) < 10 * eps)) = 0 ; % Prise en compte des erreurs de calcul
warning off
Lambda = -A./B  ;
warning on ;
% ---> Segment PiPi+1
A = dC1(:,1) * C2(1:end-1,2)' -  dC1(:,2) * C2(1:end-1,1)' + ...
    (C1(1:end-1,1) .* C1(2:end,2) - C1(1:end-1,2) .* C1(2:end,1)) * ones(1,size(dC2,1)) ;
A(find(abs(A) < 10 * eps)) = 0 ; % Prise en compte des erreurs de calcul
warning off
Kappa = -A./B ;
warning on

% c) Search the true intersections (no parallelism)
[I,J] = find((Lambda < 1)&(Lambda > 0)&(Lambda ~= NaN)&(abs(Lambda) ~= Inf) & ...
    (Kappa < 1)&(Kappa > 0)&(Kappa ~= NaN)&(abs(Kappa) ~= Inf));
%
if size(I,1) > size(I,2)
    Lambda = Lambda' ;
    tmpI = I' ;
    tmpJ = J' ;
    I = tmpJ ;
    J = tmpI ;
end

% Intersection coordinates
%
if ~isempty(I) ;
    Pts = C1(I,1:2) + repmat(Lambda(sub2ind(size(Lambda),I,J))',1,2) .* dC1(I,1:2) ;
    Ou = [I,J] ;
    if nargin == 1 ;
        % Only one vurve : spcial case for ou
        liste = find(abs(diff(Ou,1,2)) ~= 1) ;
        if ~isempty(liste) ;
            Pts = Pts(liste,:) ; Ou = Ou(liste,:) ;
        else
            Pts = [] ; Ou = [] ;
        end
    end
else
    Pts = [] ; Ou = [] ;
end

