% function [Pt2,dist] = projection_point_droite(Pts,Line) ;
%
% Projection of points on a line
%
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 24/06/2024

function [Pt2,dist] = projection_point_droite(Pt,Droite) ;

if ~isfield(Droite,'V_dir') & ~isfield(Droite,'pts') & isequal(size(Droite),[2 3])
    D = Droite ;
    Droite.V_dir = D(2,:) ;
    Droite.pts = D(1,:) ;
end

Droite.V_dir = Droite.V_dir / norm(Droite.V_dir) ;
Pt2 = repmat(Droite.pts,size(Pt,1),1) + ...
    dot(repmat(Droite.V_dir,size(Pt,1),1),Pt-repmat(Droite.pts,size(Pt,1),1),2)...
    *Droite.V_dir ;
dist = Norm2(Pt-Pt2) ;

% fini