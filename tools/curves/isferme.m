% function isclosed = isferme(Xcnt, Ycnt);
%
% Check if a set of contours are closed contours.
%
% Author: Antoine Serrurier
% Date (adaptation): 24/06/2024
%

function isclosed = isferme(Xcnt, Ycnt)

isclosed = NaN;

seuil = 10^(-10);

[nb_zones, ind_deb, ind_fin] = nb_contig_zones(Xcnt);

for i_zone = 1:nb_zones
	
	% Current sub-contour
	Xcnt_cour = Xcnt(ind_deb(i_zone):ind_fin(i_zone));
	Ycnt_cour = Ycnt(ind_deb(i_zone):ind_fin(i_zone));
	
	% Distance between extremeties
	dextr = sqrt((Xcnt_cour(end) - Xcnt_cour(1))^2 + (Ycnt_cour(end) - Ycnt_cour(1))^2);
	
	% Output
	isclosed(i_zone) = (dextr <= seuil);
	
end  % for i_zone = 1:nb_zones


















