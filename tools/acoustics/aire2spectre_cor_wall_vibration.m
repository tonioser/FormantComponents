function H_eq = aire2spectre_cor_wall_vibration(ORAL, NASAL, nbfreq, Fmax, Fmin, F_form, Z_form);

% function H = aire2spectre_cor_wall_vibration(ORAL, NASAL, nbFreq, Fmax, Fmin, F, Z);
% 
% Returns the acoustic transfer function given the oral and nasal area
% function in input where the poles (=formants) and zeros provided in input
% have been removed.
% 
% See aire2spectre_wall_vibration for further details on the inputs/outputs
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 26/06/2024

% Transfer function
H_eq = aire2spectre_wall_vibration(ORAL, NASAL, nbfreq, Fmax, Fmin);

% Frequencies
f = linspace(Fmin, Fmax, nbfreq); 
if isreal(f) SI = complex(0, 2 * pi * f); else SI = 2 * pi * j * f; end

% Remove the poles
NF = size(F_form, 1);
for I = 1 : NF
	SI1 = complex(F_form(I, 4) * pi, F_form(I, 2) * 2 * pi);
	H_eq = H_eq .* ((SI - SI1) .* (SI-conj(SI1))) ./ (SI1 .* conj(SI1));
end

% Remove the zeros
if exist('Z_form')
	NZ = size(Z_form, 1);
	for I = 1 : NZ
		SI1 = complex(Z_form(I, 4) * pi, Z_form(I, 2) * 2 * pi);
		H_eq = H_eq .* (SI1 .* conj(SI1)) ./ ((SI - SI1) .* (SI-conj(SI1))) ;
	end
end

return
