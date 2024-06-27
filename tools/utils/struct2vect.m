function vect = struct2vect(STRUCT, ind_struct, var, ind_var)

% function vect = struct2vect(STRUCT, ind_struct, var, ind_var);
%
% Tranforms data in structures into a vector
%
% Inputs
%     STRUCT(struct)       : Input structures
%     ind_struct(nb_struct): Selected structures from the input structures
%     var(string)          : Name of the selected variable in the structures
%     ind_var(nb_var)      : Indices of the selected indices in the variable 
%
% outputs
%     vect: Output vector
%
% Example:
%     vect = struct2vect(CONF, [1,2,10], 'error', [2,3]);
%     => vect == [CONF(1).erreur(2), CONF(1).error(3), CONF(2).error(2),...
%                 CONF(2).erreur(3), CONF(10).error(2), CONF(10).error(3)]
% 
% Author: Antoine Serrurier
% Date (adaptation): 26/06/2024


% Initialisation
vect = [];

% Loop on the selected structures
for i_struct = ind_struct
	% Loop on the indices in thz variable
	for i_var = ind_var
		vect = [vect, eval(['STRUCT(i_struct).', var, '(i_var)'])];
	end  % for i_var = ind_var
end  % for i_struct = ind_struct











