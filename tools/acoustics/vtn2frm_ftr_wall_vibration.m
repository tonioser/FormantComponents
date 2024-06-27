function [H_eq, F_form, Z_form] = vtn2frm_ftr_wall_vibration(ORAL, NASAL, nbfreq, Fmax, Fmin);
%
% [H, F, Z] = vtn2frm_ftr_wall_vibration(ORAL, NASAL, nbFreq, Fmax, Fmin);
% 
% Returns the acoustic transfer function and the values of the formants
% and zeros given the area functions of the oral and nasal tracts in input
% 
% Inputs
%   ORAL      : Structure containing the parameters of the oral area function
%       ORAL.A(nbOralTubes)         : Areas of the tubes in cm2 (glottis -> lips)
%       ORAL.L(nbOralTubes)	        : Lengths of the tubes in cm (glottis -> lips)
%
%   NASAL     : Structure containing the parameters of the nasal area function
%       NASAL.A_mid(nbNasophaTubes) : Areas of the tubes of the nasopharyngeal port in cm2 
%       NASAL.L_mid(nbNasophaTubes) : Lengths of the tubes of the nasopharyngeal port in cm
%       NASAL.ind_connect_mid       : Index in the oral area function of the connection of the nasopharyngeal port 
%
%       NASAL.A_left(nbNasLeft)     : Areas of the tubes of the left nasal cavity in cm2 
%       NASAL.L_left(nbNasLeft)     : Lengths of the tubes of the left nasal cavity in cm
%       NASAL.A_right(nbNasRight)   : Areas of the tubes of the right nasal cavity in cm2 
%       NASAL.L_right(nbNasRight)   : Lengths of the tubes of the right nasal cavity in cm
%
%       NASAL.A_sinleft(nbSinLeft)  : Areas of the tubes of the left maxillary sinus in cm2 
%       NASAL.L_sinleft(nbSinLeft)  : Lengths of the tubes of the left maxillary sinus in cm
%       NASAL.A_sinright(nbSinRight): Areas of the tubes of the right maxillary sinus in cm2 
%       NASAL.L_sinright(nbSinRight): Lengths of the tubes of the right maxillary sinus in cm
%       NASAL.ind_connect_sinleft   : Index in the left nasal cavity area function of the connection of the left maxillary sinus 
%       NASAL.ind_connect_sinright  : Index in the left right cavity area function of the connection of the right maxillary sinus
%   nbFreq(1) : Number of frequency points calculated (optional, [default = 256])
%   Fmax(1)   : Maximal frequency calculated in Hz (optional, [default = 5000])
%   Fmin(1)   : Minimal frequency calculated in Hz (optional, [default = Fmax / nbfreq])
%
% Outputs
%     H(nbFreq) : Transfer function (amplitude plot: '20*log10(abs(H))')
%     F(nbForm) : Structure containing formant values (poles of the transfer function)
%     Z(nbForm) : Structure containing zero values (zeros of the transfer function)
%
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 25/06/2024


%===================================================================
% Inputs and constants

% Input arguments
if nargin < 3
    nbfreq = 256; Fmax = 5000; Fmin = Fmax / nbfreq;
end  % if nargin < 3

% Constant
seuil2 = 10;

%  Vector of frequencies
f = linspace(Fmin, Fmax, nbfreq);

%===================================================================
% Acoustic transfer function, taking into account wall vibration

H_eq = aire2spectre_wall_vibration(ORAL, NASAL, nbfreq, Fmax, Fmin);

%===================================================================
% Formant (poles) and zero detection

%------------------------------------------------------------------------
% Peaks of the transfer function according to peakpick (rough formants)

[maxi, ind_maxi, mini, ind_mini] = peakpick(20*log10(abs(H_eq)));
frq_min = f(ind_mini);
frq_max = f(ind_maxi);
FEST = frq_max;
IFEST = length(FEST);

%------------------------------------------------------------------------
% Refinement according to nraph

FE = 300; BNPE = 50; ITERMX = 100; FMAX = 5000; NF = 0; F_form = []; Z_form = [];
for II = 1 : IFEST
    FE = FEST(II);
    % Look for the poles of H_eq by nraph
    [F, BP] = nraph_wall_vibration(FE, BNPE, ITERMX, FMAX, ORAL, NASAL, F_form);
    NF = NF + 1; F_form(NF, 2) = F; F_form(NF, 4) = BP;
    % Discard if we find the same pole
    if NF > 1
        if (abs(F_form(NF, 2) - F_form(NF-1, 2))) < seuil2
            F_form(NF, :) = []; NF = NF - 1;
        end
    end % if NF > 1
end

%------------------------------------------------------------------------
% Recalculation of the acoustic transfer function by removing the poles of the already detected formants
% This might reveal poles that were "hidden" by the dominent poles

H_eq_cor = aire2spectre_cor_wall_vibration(ORAL, NASAL, nbfreq, Fmax, Fmin, F_form);

%------------------------------------------------------------------------
% Peaks of the second transfer function according to peakpick (rough formants)

[maxi, ind_maxi, mini, ind_mini] = peakpick(20*log10(abs(H_eq_cor)));
frq_min = f(ind_mini); frq_max = f(ind_maxi);
FEST = frq_max; IFEST = length(FEST);

%------------------------------------------------------------------------
% Refinement according to nraph

for II = 1 : IFEST
    FE = FEST(II);
    % Look for the poles of H_eq by nraph
    [F, BP] = nraph_wall_vibration(FE, BNPE, ITERMX, FMAX, ORAL, NASAL, F_form);
    NF = NF + 1; F_form(NF, 2) = F; F_form(NF, 4) = BP;
    % Discard if we find the same pole
    if NF > 1
        if (abs(F_form(NF, 2) - F_form(NF-1, 2))) < seuil2
            F_form(NF, :) = []; NF = NF - 1;
        end
    end % if NF > 1
end

%------------------------------------------------------------------------
% Recalculation, again, of the acoustic transfer function by removing the poles of the already detected formants
% This may reveal zeros

H_eq_cor = aire2spectre_cor_wall_vibration(ORAL, NASAL, nbfreq, Fmax, Fmin, F_form);

%------------------------------------------------------------------------
% Mins of the third transfer function according to peakpick

[maxi, ind_maxi, mini, ind_mini] = peakpick(20*log10(abs(H_eq_cor)));
frq_min = f(ind_mini); frq_max = f(ind_maxi);
ZEST = frq_min; IZEST = length(ZEST);

%------------------------------------------------------------------------
% Refinement of the zeros according to nraph

FE = 300; BNPE = 50; NZ = 0; 
for II = 1 : IZEST
	FE = ZEST(II);
    % Look for the poles of H_eq by nraph
  [Z, BZ] = nraph_wall_vibration(FE, BNPE, ITERMX, FMAX, ORAL, NASAL, F_form, Z_form);
	NZ = NZ + 1; Z_form(NZ, 2) = Z; Z_form(NZ, 4) = BZ;
    % Discard if we find the same pole
	if NZ > 1
		if (abs(Z_form(NZ, 2) - Z_form(NZ-1, 2))) < seuil2
			Z_form(NZ, :) = []; NZ = NZ - 1;
		end
	end % if NZ > 1
end

%===================================================================
% Poles (=formants) and zeros cleaning
%
% Removing duplicate poles/zeros
% Sorting
% Removing poles/zeros beyond Fmax
% Removing poles/zeros with negative frequencies

% Poles
if ~isempty(F_form)
	% Sort
	[bid, ind_F_ord] = sort(F_form(:,2));
	F_form = F_form(ind_F_ord,:);
	% Remove poles below 0 and beyond Fmax
	F_form = F_form(find((F_form(:,2) <= Fmax) & (F_form(:,2) >= 0)),:);
	% Remove duplicate poles
	ind_F_diffOK = find(diff(F_form(:,2)) >= seuil2);
	if ~isempty(ind_F_diffOK)
		ind_F_Nredondant = [ind_F_diffOK; ind_F_diffOK(end) + 1];
	elseif isempty(ind_F_diffOK) & ~isempty(F_form) % if ~isempty(ind_F_diffOK)
		ind_F_Nredondant = 1;
	else  % if ~isempty(ind_F_diffOK)
		ind_F_Nredondant = [];
	end  % if ~isempty(ind_F_diffOK)
	F_form = F_form(ind_F_Nredondant,:);
end  % if ~isempty(F_form)
if isempty(F_form); F_form = NaN * ones(1,4); end

% Zeros
if ~isempty(Z_form)
	% Sort
	[bid, ind_Z_ord] = sort(Z_form(:,2));
	Z_form = Z_form(ind_Z_ord,:);
	% Remove zeros below 0 and beyond Fmax
	Z_form = Z_form(find((Z_form(:,2) <= Fmax) & (Z_form(:,2) >= 0)),:);
	% Remove duplicate poles
	ind_Z_diffOK = find(diff(Z_form(:,2)) >= seuil2);
	if ~isempty(ind_Z_diffOK)
		ind_Z_Nredondant = [ind_Z_diffOK; ind_Z_diffOK(end) + 1];
	elseif isempty(ind_Z_diffOK) & ~isempty(Z_form) % if ~isempty(ind_Z_diffOK)
		ind_Z_Nredondant = 1;
	else  % if ~isempty(ind_Z_diffOK)
		ind_Z_Nredondant = [];
	end  % if ~isempty(ind_Z_diffOK)
	Z_form = Z_form(ind_Z_Nredondant,:);
end  % if ~isempty(Z_form)
if isempty(Z_form); Z_form = NaN * ones(1,4); end

%===================================================================
% Amplitude of the transfer function, in dB, for the formants and zeros

for i_pole = 1:size(F_form,1)
	F_form(i_pole,3) = 20*log10(abs(aire2spectre_wall_vibration(ORAL, NASAL, 1, F_form(i_pole,2), F_form(i_pole,2))));
end  % for i_pole = 1:size(F_form,1)

for i_zero = 1:size(Z_form,1)
	Z_form(i_zero,3) = 20*log10(abs(aire2spectre_wall_vibration(ORAL, NASAL, 1, Z_form(i_zero,2), Z_form(i_zero,2))));
end  % for i_pole = 1:size(F_form,1)

return







