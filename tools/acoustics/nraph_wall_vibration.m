% function [F, BP] = nraph_wall_vibration(FE, BNPE, ITERMX, FMAX, ORAL, NASAL, F_form, Z_form);
% 
% Returns the poles/zeros frequencies and bandwidth in an acoustic transfer
% function given the oral and nasal area functions and an estimation of the
% locations of the poles/zeros.
% 
% See aire2spectre_wall_vibration for further details about the
% input/outputs.
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 26/06/2024



function [F, BP, ITER, FRQ, HHH] = nraph_wall_vibration(FE, BNPE, ITERMX, FMAX, ORAL, NASAL, F_form, Z_form);

% Initialisations
DELTAS = complex(30, 30); SEUIL = 0.3; ITER = 1;
SI = complex(-BNPE*pi, FE*2*pi);

F = NaN; BP = NaN;

% Depending if Z_form is provided we look for the poles or the zeros

% Iteration until ITERMX
while ITER < ITERMX
	% First point
	FIcx = -j * SI / (2*pi); 
	if ~exist('Z_form')
		Q = 1 / aire2spectre_cor_wall_vibration(ORAL, NASAL, 1, FIcx, FIcx, F_form);
	else
		Q = aire2spectre_cor_wall_vibration(ORAL, NASAL, 1, FIcx, FIcx, F_form, Z_form);
	end

	% Second point
	SIP = SI + DELTAS; FIPcx = -j * SIP / (2*pi);
	if ~exist('Z_form')
		QP = 1 / aire2spectre_cor_wall_vibration(ORAL, NASAL, 1, FIPcx, FIPcx, F_form);
	else
		QP = aire2spectre_cor_wall_vibration(ORAL, NASAL, 1, FIPcx, FIPcx, F_form, Z_form);
	end
	% intersection point with 0
	Q1D = (QP - Q) / DELTAS; SISU = SI - (Q / Q1D);
	
	% If found...
	if(abs(SISU-SI) < SEUIL) F = imag(SISU) / (2*pi); BP = - real(SISU)/pi; return
	% ... otherwise
	else SI = SISU; ITER = ITER + 1;	end
end  % while ITER < ITERMX	


return
