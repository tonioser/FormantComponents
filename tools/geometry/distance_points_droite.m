% dist = distance_points_droite(line,pts) 
%
% Fonction de calcul de distance entre une droite et des points
% Returns the distances between points and a line
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 24/06/2024

function [dist,dpto] = distance_points_droite(D,P) ;

% Input arguments
if ~isstruct(D) ;
    DD.V_dir = D(2,:) ;
    DD.pts = D(1,:) ;
    D = DD ;
% else
%     Donnees = {varargin{:}} ;
end

if (length(D) ~= size(P,2))&(length(D) ~= 1) ;
    % ---> Error : no calculation possible
    error('Wrong input format') ;
end

% 2. Variable formatting
dim = size(P,2) ; % ---> Dimension
N = size(P,1) ;   % ---> Number of points
if length(D) == 1 ;
    % Matrice of the origin points and direction vectors
    Pto = ones(N,1) * D.pts ;
    MVd = ones(N,1) * D.V_dir ;
else
    % Matrice of the origin points and direction vectors
    Pto = reshape([D(:).pts],Dim,N)' ;
    MVd = reshape([D(:).V_dir],Dim,N)' ;
end
% All the direction vectors
MVd = norme_vecteur(MVd) ;

% 3. Distances
dpto = dot(MVd,P-Pto,2) ; % distance of the projection of the points on the line
dist = sqrt(Norm2(P-Pto).^2 - dpto.^2) ;
