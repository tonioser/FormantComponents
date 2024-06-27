function [DistSag, Long, GRD_out] = crossarea_2_tract_2D(GRD)

%
% function [DistSag, Long, GRD_out] = crossarea_2_tract_2D(GRD);
%
% Returns the sagittal distance given the intersections between the vocal
% tract outlines and a grid
%
% Inputs
%     GRD(nb_grd).X(1,2) and Y(1,2)                 : Coordinates of the intersection points between the tube the and grid lines
%     GRD(nb_grd).X2D_grd_Int(1) and Y2D_grd_Int(1) : Interior point for erach grid line 
%     GRD(nb_grd).X2D_grd_Ext(1) and Y2D_grd_Ext(1) : Exterior point for erach grid line 
%
% outputs
%     DistSag(1,nb_grd-1) : Sagittal distances for the (nb_grd-1) tubes of the sagittal function
%     Long(nb_grd-1)    : Lengths for the (nb_grd-1) tubes of the sagittal function
%     GRD_out(nb_grd)   : Input structure completed with sagittal distances, gravity centres and correction factor for each grid line
%
% Author: Antoine Serrurier
% Date: 13/07/2009

% Number of grid lines
nb_grd = length(GRD);

%---------------------------------------------------------------------
% Distances and gravity centres of the grid lines
% Expansion of the gravity centres in 3D

for i_grd = 1:nb_grd
	
	% Initialisation of the distances and gravity centres
	%-----------------------------------------------
	% Distance
	GRD(i_grd).D = 0;
	% Gravity centres
	GRD(i_grd).Xc = 0;
	GRD(i_grd).Yc = 0;
	
	% Distances and gravity centres
	%---------------------------------------
	if ~all(isnan(GRD(i_grd).X))
		
		% Even number of points
		if mod(length(GRD(i_grd).X), 2) == 1
			% If only one point, duplicate it
			if length(GRD(i_grd).X) == 1
				GRD(i_grd).X = [GRD(i_grd).X, GRD(i_grd).X];
				GRD(i_grd).Y = [GRD(i_grd).Y, GRD(i_grd).Y];
				% errordlg(['WARNING crossarea_2_tract_2D: the single point of the grid ', num2str(i_grd), ' has been duplicated!'])
			else  % if length(GRD(i_grd).X) == 1
				error(['ERROR crossarea_2_tract_2D: the number of points of the grid ', num2str(i_grd), ' is odd!'])
				errordlg(['ERROR crossarea_2_tract_2D: the number of points of the grid ', num2str(i_grd), ' is odd!'])
			end  % if length(GRD(i_grd).X) == 1
		end  % if mod(length(GRD(i_grd).X), 2) == 1
		
		% Point sorting from lower to upper
		[ii, jj] = max(dist([GRD(i_grd).X; GRD(i_grd).Y]', [GRD(i_grd).X; GRD(i_grd).Y]));
		ind_pt_org = jj(1);
		[bid, ind_pt_ord] = sort(dist([GRD(i_grd).X(ind_pt_org); GRD(i_grd).Y(ind_pt_org)]', [GRD(i_grd).X; GRD(i_grd).Y]));
		
		% Distance between consecutive points
		DDtmp = dist([GRD(i_grd).X(ind_pt_ord(1:end-1)); GRD(i_grd).Y(ind_pt_ord(1:end-1))]',...
		          [GRD(i_grd).X(ind_pt_ord(2:end)); GRD(i_grd).Y(ind_pt_ord(2:end))]);
		DD = DDtmp(1:(size(DDtmp,1)+1):end);
		
		% Positive distances
		D_plus = sum(DD(1:2:end));
		
		% Negative distances
		D_moins = sum(DD(2:2:end));
		
		% Gravity centres of the positive regions
		Xc_plus = mean([GRD(i_grd).X(ind_pt_ord(1:2:end-1));GRD(i_grd).X(ind_pt_ord(2:2:end))]);
		Yc_plus = mean([GRD(i_grd).Y(ind_pt_ord(1:2:end-1));GRD(i_grd).Y(ind_pt_ord(2:2:end))]);
		
		% Final weighted gravity centre
		Xc = sum(DD(1:2:end) .* Xc_plus) ./ sum(DD(1:2:end));
		Yc = sum(DD(1:2:end) .* Yc_plus) ./ sum(DD(1:2:end));
		
		%**********************************************
		while 0
		figure
		FIG
		plot(GRD(i_grd).X, GRD(i_grd).Y, 'r*')
		text(GRD(i_grd).X(ind_pt_ord), GRD(i_grd).Y(ind_pt_ord), num2str((1:length(GRD(i_grd).X))'))
		plot(Xc_plus, Yc_plus, 'b*')
		plot(Xc, Yc, 'k*')
		end  % while 0
		%**********************************************
		
		% Save
		GRD(i_grd).D = D_plus;
		GRD(i_grd).Xc = Xc;
		GRD(i_grd).Yc = Yc;
		
	else  % if ~all(isnan(GRD(i_grd).PLUS.X))
		GRD(i_grd).D = NaN;
		GRD(i_grd).Xc = NaN;
		GRD(i_grd).Yc = NaN;
	end  % if ~all(isnan(GRD(i_grd).PLUS.X))
	
end  % for i_grd = 1:nb_grd

%------------------------------------------------------------------
% Lengths

% Simple calculation
Long = NaN * ones(1,nb_grd-1);
for i_grd = 1:nb_grd-1
	if (~isnan(GRD(i_grd).Xc)) & (~isnan(GRD(i_grd+1).Xc)) & (GRD(i_grd).D ~= 0) & (GRD(i_grd+1).D ~= 0)
		Long(i_grd) = dist([GRD(i_grd).Xc, GRD(i_grd).Yc], [GRD(i_grd+1).Xc, GRD(i_grd+1).Yc]');
	end  % if (~isnan(GRD(i_grd).Xc3D)) & (~isnan(GRD(i_grd+1).Xc3D)) & (GRD(i_grd).A ~= 0) & (GRD(i_grd+1).A ~= 0)
end  % for i_grd = 1:nb_grd-1

% Treat the case where the area is 0
% If the area is 0, no gravity centre => extrapolated from the neighbours 
% Series areas where all areas are equal to 0
aire = struct2vect(GRD,1:nb_grd,'D',1);
ind_zeros_a_traiter = find(aire == 0);
% Loop on the number of regions with area equal to 0
while ~isempty(ind_zeros_a_traiter)
	ind_zeros = [ind_zeros_a_traiter(1)];
	% Look for several consecutive
	diff_ind = diff(ind_zeros_a_traiter);
	if isempty(diff_ind); diff_ind = 1000; end;
	while (diff_ind(1) == 1)
		ind_zeros = [ind_zeros, ind_zeros(end) + 1];
		diff_ind = diff_ind(2:end);
		if isempty(diff_ind); diff_ind = 1000; end;
	end  % while (diff_ind(1) == 1)
	
	% Length between the two nearest non-null CG
    % If nothing at the end end, we set a length equal to the previous length
	if length(GRD) >= ind_zeros(end)+1
		if ~isnan(GRD(ind_zeros(end)+1).Xc)
			L = dist([GRD(ind_zeros(1)-1).Xc, GRD(ind_zeros(1)-1).Yc],...
			         [GRD(ind_zeros(end)+1).Xc, GRD(ind_zeros(end)+1).Yc]');
		else  % if ~isnan(GRD(ind_zeros(end)+1).ORAL.Xc3D)
			L = Long((ind_zeros(1)-2)) * (length(ind_zeros) + 1);
		end  % if ~isnan(GRD(ind_zeros(end)+1).Xc3D)
	else
		L = Long((ind_zeros(1)-2)) * (length(ind_zeros) + 1);
	end  % if length(GRD) >= ind_zeros(end)+1
		
    % Split the length by the number of considered tubes
	for i_grd = ind_zeros
		Long(i_grd-1) = L / (length(ind_zeros) + 1);
		Long(i_grd) = L / (length(ind_zeros) + 1);
	end  % for i_grd = ind_zeros
	
	% Updates cases remaining to treat
	[bid, ind_non_traites] = vectcmp(ind_zeros_a_traiter', ind_zeros');
	ind_zeros_a_traiter = ind_zeros_a_traiter(ind_non_traites);
end  % while ~isempty(ind_zeros_a_traiter)

%------------------------------------------------------------------
% Areas

% Principe: pour corriger une éventuelle coupe par un plan de biais par
% rapport à la perpendiculaire au tuyau du conduit vocal qui conduirait
% à surévaluer l'aire calculée dans ce plan (coupe en baiais d'un soucisson
% => aire d'une ellipse au lieu de l'aire d'un cercle), on corrige
% l'aire de chaque plan de grille.
%
% Calculs:
%
% 1. Pour estimer la forme du tuyau en un point de mesure, on l'approxime
%    par un cercle qui passe par le centre de gravité considéré et celui
%    de ses 2 grilles voisines, dans le plan formé par ces 3 centres
%    de gravité
% 2. la droite passant par le centre du cercle et le centre de gravité
%    considéré constitue l'intersection entre le plan de coupe idéal
%    et le plan des 3 centres de gravité (qui lui est perpendiculaire)
% 3. On corrige alors l'aire par le cosinus de l'angle 3D entre le vecteur
%    normal au plan de coupe et le vecteur normal au plan de coupe idéal.

% Expand to 3D to re-use a code originally written for the 3D case 
clear GRD3D Plan_grd_cmicp
for i_grd = 1:nb_grd
	GRD3D(i_grd).Xc3D = GRD(i_grd).Xc;
	GRD3D(i_grd).Yc3D = 0;
	GRD3D(i_grd).Zc3D = GRD(i_grd).Yc;
	
	Plan_grd_cmicp(:,1,i_grd) = [GRD(i_grd).X2D_grd_Int, 0, GRD(i_grd).Y2D_grd_Int];
	Plan_grd_cmicp(:,2,i_grd) = [GRD(i_grd).X2D_grd_Ext, 0, GRD(i_grd).Y2D_grd_Ext];
	Plan_grd_cmicp(:,3,i_grd) = [GRD(i_grd).X2D_grd_Int, -1, GRD(i_grd).Y2D_grd_Int];
end  % for i_grd = 1:nb_grd

% Correction factor 
for i_grd = 1:nb_grd
	
	% the 3 gravity centres
    % If one is missing, defined as the continuation of the others 
	if (i_grd == 1) || (isnan((GRD3D(i_grd-1).Xc3D)) & ~(i_grd == nb_grd))
		CG2 = [GRD3D(i_grd).Xc3D, GRD3D(i_grd).Yc3D, GRD3D(i_grd).Zc3D];
		CG3 = [GRD3D(i_grd+1).Xc3D, GRD3D(i_grd+1).Yc3D, GRD3D(i_grd+1).Zc3D];
		CG1 = CG2 + diff([CG3; CG2]);
	elseif (i_grd == nb_grd) || isnan((GRD3D(i_grd+1).Xc3D))  % if (i_grd == 1) || isnan((GRD3D(i_grd-1).Xc3D))
		CG1 = [GRD3D(i_grd-1).Xc3D, GRD3D(i_grd-1).Yc3D, GRD3D(i_grd-1).Zc3D];
		CG2 = [GRD3D(i_grd).Xc3D, GRD3D(i_grd).Yc3D, GRD3D(i_grd).Zc3D];
		CG3 = CG2 + diff([CG1; CG2]);
	else  % if (i_grd == 1) || isnan((GRD3D(i_grd-1).Xc3D))
		CG1 = [GRD3D(i_grd-1).Xc3D, GRD3D(i_grd-1).Yc3D, GRD3D(i_grd-1).Zc3D];
		CG2 = [GRD3D(i_grd).Xc3D, GRD3D(i_grd).Yc3D, GRD3D(i_grd).Zc3D];
		CG3 = [GRD3D(i_grd+1).Xc3D, GRD3D(i_grd+1).Yc3D, GRD3D(i_grd+1).Zc3D];
	end  % if (i_grd == 1) || isnan((GRD3D(i_grd-1).Xc3D))
	
	%*********************************
	while 0
	figure
	FIG
	PLOT3([CG1(1), CG2(1), CG3(1)], [CG1(2), CG2(2), CG3(2)], [CG1(3), CG2(3), CG3(3)], '*-')
	text([CG1(1), CG2(1), CG3(1)], [CG1(2), CG2(2), CG3(2)], [CG1(3), CG2(3), CG3(3)], str2mat('CG1','CG2','CG3'))
	end  % while 0
	%*********************************
	
	% Si l'un des 3 centres de gravité est nul, cela signifie que l'aire de la section
	% correspondante est nulle, que l'aire des tubes de part et dautre de cette section
	% est aussi nulle et donc qu'il ne sert à rien de calculer cosalp (il n'est d'ailleurs
	% pas possible d'estimer la direction du tuyau)

    % If one of the 3 CG is null, it means that the area of the
    % cross-section is equal to 0 and the area of one of the tube around
    % the grid in 0 and there is no need to calculate the correction factor
	if all(CG1 == [0 0 0]) || all(CG2 == [0 0 0]) || all(CG3 == [0 0 0])
		cosalp(i_grd) = 0;
	else  % if all(CG1 == [0 0 0]) || all(CG2 == [0 0 0]) || all(CG3 == [0 0 0])
		cosalp(i_grd) = corrige_crossarea_oblique(CG1, CG2, CG3, squeeze(Plan_grd_cmicp(:,:,i_grd)));
	end  % if all(CG1 == [0 0 0]) || all(CG2 == [0 0 0]) || all(CG3 == [0 0 0])
	
	GRD(i_grd).cosalp = cosalp(i_grd);
	
end  % for i_grd = 1:nb_grd

% Distances as the weight of the 2 corrected distances
DistSag = ((aire(1:end-1) .* cosalp(1:end-1)) + (aire(2:end) .* cosalp(2:end))) / 2;

% Output
GRD_out = GRD;
	
return






