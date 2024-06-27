function cosalp = corrige_crossarea_oblique(CG1, CG2, CG3, Plan_cmicp)

% function cosalp = corrige_crossarea_oblique(CG1, CG2, CG3, Plane);
%
% Returns the correction factor to apply to the surface of the
% cross-section of a tube cut by a plane not exactly perpendicular to the
% tube direction.
% The cross-section surface needs to be multiplied by this factor to
% estimate the cross-section surface if the plane were exactly perpendiclar
% to the direction of the tube.
%
% The direction of the tube is estimated by a circle passing by the centre
% of cross section the centres of gravity of the neighbouring
% cross-sections.
% 
%
% Inputs
%     CG1(1,iXYZ) and CG3(1,iXYZ) : Coordinates of the gravity centres of the neighbouring cross-sections 
%     CG2(1,iXYZ)                 : Coordinates of the gravity centre of the current cross-section 
%     Plane(iXYZ,3)               : Coordinates of three points defining the cutting plane 
%
% outputs
%     cosalp(1) : Correction factor (<1)
%
% Author: Antoine Serrurier
% Date (adapted): 26/06/2024

% Explication mathématique:
%
% Considérons le plan formé par les 3 centres de gravités. On travaille
% en 2D dans ce plan.
% On cherche à caractériser la droite (resp. le plan) qui passe par
% le point CG2 et le centre du cercle passant par CG1, CG2 et CG3.
% - Soit lambda l'angle entre les 2 segments [CG1 CG2] et [CG2 CG3];
%   c'est aussi l'angle entre leur médiatrices respectives.
%   Cet angle est connu.
% - Soit alpha l'angle entre la droite recherchée et le segment
%   [CG1 CG2]. Cet angle est inconnu.
% - Soit beta l'angle entre la droite recherchée et le segment
%   [CG2 CG3]. Cet angle est inconnu.
% - Soit d1 et d2 les demi-longueurs respectices des segments
%   [CG1 CG2] et [CG2 CG3]
% On peut alors montrer dans le plan les relations suivantes:
%  1. alpha + beta + lambda = pi
%  2. d1 cos(beta) = d2 cos (alpha)
% On en déduit alors:
% tan(alpha) = (d2 + d1 cos(lambda)) / (d1 sin(lambda))
% alpha suffit à caractériser complètement la droite (resp. le plan)
% recherchée en 2D (resp. en 3D).
% On généralise facilement en 3D.
%
% Le facteur de correction est alors le cosinus de l'angle entre
% la normale du plan idéal (cad celui dont on vient de calculer
% la position en considérant que le tuyau était de direction
% circulaire) et la normale du plan de coupe réel utilisé.

% Gravity centres, defining a plane
CG = [CG1; CG2; CG3]';

% Half-lengths of each segment
d1 = dist(CG1, CG2') / 2;
d2 = dist(CG2, CG3') / 2;

% Angle between the segments
seg1 = CG2 - CG1;
seg2 = CG3 - CG2;
lambda = real(acos(dot(seg1,seg2) / (norm(seg1) * norm(seg2))));

% Angle between the normal to the ideal plane and 1
alp = atan2(d2 + d1 * cos(lambda), d1 * sin(lambda));

% Normal to the ideal plane = segment 1 rotated of pi / 2 - alp in the plane (CG1 CG2 CG3)
% Check if (CG1 CG2 CG3) form a plane (or, equivalently, if alp ~= pi / 2, which means no rotation)
if norm(cross(seg1, seg2)) ~= 0
	CG12_rot = rot3d_plan([CG1 ; CG2], CG, pi / 2 - alp, CG1);
else  % if norm(cross(seg1, seg2)) ~= 0
	CG12_rot = [CG1 ; CG2];
end  % if norm(cross(seg1, seg2)) ~= 0

% Normal of the ideal plance
norm_ideal = diff(CG12_rot);
norm_ideal = norm_ideal / norm(norm_ideal);

% Normal of the current plane
norm_cpe = cross(Plan_cmicp(:,3) - Plan_cmicp(:,1), Plan_cmicp(:,2) - Plan_cmicp(:,1))';
norm_cpe = norm_cpe / norm(norm_cpe);

% Cosine of the 3D angle between these 2 normals
cosalp = abs(dot(norm_cpe, norm_ideal));

return





