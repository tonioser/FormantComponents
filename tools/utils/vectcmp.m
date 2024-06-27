% function [ind_equal, ind_diff] = vectcmp(vect1, vect2 {, seuil_distance});
%
% Returns the indices of vectir vect1 having similar and different values
% than vector vect2.
%
% A threshold can be set under which 2 values are considered as equal.
%
% Inputs
%   vect1(n,1) : Input vector 1
%   vect2(p,1) : Input vector 2
% 
% Outputs
%   ind_equal(1,q) : Indices in vect1 of values present in vect2
%   ind_diff(1,q)  : Indices in vect1 of values not present in vect2
%
% Author: Antoine Serrurier
% Date: 27/05/2005

function [ind_equal, ind_diff] = vectcmp(vect1, vect2, seuil_distance)

%----------------------------------------------------------
% Trivial case when vect1 is a single value

if (length(vect1) == 1)

    if nargin == 3
        test_bool = ~isempty(find(abs(vect1 - vect2) <= seuil_distance));
    else  % if nargin == 3
        test_bool = ~isempty(find(vect1 == vect2));
    end  % if nargin == 3

    if test_bool
        ind_equal = 1; ind_diff = [];
    else  % if ~isempty(find(vect1 == vect2))
        ind_equal = []; ind_diff = 1;
    end  % if ~isempty(find(vect1 == vect2))

%----------------------------------------------------------
% Trivial case when vect2 is a single value

elseif (length(vect2) == 1)  % if (length(vect1) == 1)

    if nargin == 3
        ind_equal = find(abs(vect1 - vect2) <= seuil_distance)';
        ind_diff = find(abs(vect1 - vect2) > seuil_distance)';
    else  % if nargin == 3
        ind_equal = find(vect1 == vect2)';
        ind_diff = find(vect1 ~= vect2)';
    end  % if nargin == 3

%----------------------------------------------------------
% General case

else  % if (length(vect1) == 1)

    % Split the problem if matrices are too big => lasts longer but avoid memery issues 
    % Empirical maximal matrix size
    mat_lim_size = 10^7;

    % Split the problem if matrix size exceeded
    if (length(vect1) * length(vect2)) > mat_lim_size
        
        % Number of splits in a vector
        nb_new_vect = ceil(length(vect1) * length(vect2) / mat_lim_size) * 10;
        % Longest vector => pplitting by the longest vector
        [bid, i_vect_max] = max([length(vect1), length(vect2)]);
        
        % Longest vector = vect1
        if i_vect_max == 1

            taille_vecteurs = ceil(length(vect1) / nb_new_vect);
            nb_new_vect = ceil(length(vect1) / taille_vecteurs);
            taille_dernier_vecteur = length(vect1) - taille_vecteurs * (nb_new_vect - 1);
            for i_vect = 1:nb_new_vect-1
                ind_deb = i_vect * taille_vecteurs - taille_vecteurs + 1;
                ind_fin = ind_deb + taille_vecteurs - 1;
                VECT1(i_vect).ind = ind_deb:ind_fin;
                VECT1(i_vect).vect = vect1(ind_deb:ind_fin);
            end  % for i_vect = 1:nb_new_vect-1
            ind_deb = (nb_new_vect-1) * taille_vecteurs + 1;
            ind_fin = ind_deb + taille_dernier_vecteur - 1;
            VECT1(nb_new_vect).ind = ind_deb:ind_fin;
            VECT1(nb_new_vect).vect = vect1(ind_deb:ind_fin);

            % Call the exact same function for each subvector
            for i_vect = 1:length(VECT1)
                if nargin == 3
                    [VECT1(i_vect).ind_equal, VECT1(i_vect).ind_diff] = vectcmp(VECT1(i_vect).vect, vect2, seuil_distance);
                else  % if nargin == 3
                    [VECT1(i_vect).ind_equal, VECT1(i_vect).ind_diff] = vectcmp(VECT1(i_vect).vect, vect2);
                end  % if nargin == 3
            end  % for i_vect = 1:length(VECT2)

            % Gather outputs for all sub-vectors
            % Equal
            ind_equal = [];
            for i_vect = 1:length(VECT1)
                ind_equal = [ind_equal, (VECT1(i_vect).ind(1) - 1 + VECT1(i_vect).ind_equal)];
            end  % for i_vect = 1:length(VECT1)
            % Diff
            if ~isempty(ind_equal)
                %[bid, ind_diff] = vectcmp((1:length(vect1))', ind_equal');
                ind_tot = 1:length(vect1); ind_tot(ind_equal) = zeros(1,length(ind_equal));
                ind_diff = find(ind_tot ~= 0);
            else  % if ~isempty(ind_equal)
                ind_diff = (1:length(vect1));
            end  % if ~isempty(ind_equal)

        % Longest vector = vect2
        else  % if i_vect_max == 1

            taille_vecteurs = ceil(length(vect2) / nb_new_vect);
            nb_new_vect = ceil(length(vect2) / taille_vecteurs);
            taille_dernier_vecteur = length(vect2) - taille_vecteurs * (nb_new_vect - 1);
            for i_vect = 1:nb_new_vect-1
                ind_deb = i_vect * taille_vecteurs - taille_vecteurs + 1;
                ind_fin = ind_deb + taille_vecteurs - 1;
                VECT2(i_vect).ind = ind_deb:ind_fin;
                VECT2(i_vect).vect = vect2(ind_deb:ind_fin);
            end  % for i_vect = 1:nb_new_vect-1
            ind_deb = (nb_new_vect-1) * taille_vecteurs + 1;
            ind_fin = ind_deb + taille_dernier_vecteur - 1;
            VECT2(nb_new_vect).ind = ind_deb:ind_fin;
            VECT2(nb_new_vect).vect = vect2(ind_deb:ind_fin);

            % Call the exact same function for each subvector
            for i_vect = 1:length(VECT2)
                if nargin == 3
                    [VECT2(i_vect).ind_equal, VECT2(i_vect).ind_diff] = vectcmp(vect1, VECT2(i_vect).vect, seuil_distance);
                else  % if nargin == 3
                    [VECT2(i_vect).ind_equal, VECT2(i_vect).ind_diff] = vectcmp(vect1, VECT2(i_vect).vect);
                end  % if nargin == 3
            end  % for i_vect = 1:length(VECT2)

            % Gather outputs for all sub-vectors
            % Equal
            ind_equal_tmp = sort([VECT2.ind_equal]);
            if ~isempty(ind_equal_tmp)
                ind_equal = find([1, diff(ind_equal_tmp)] ~= 0);
                % Diff
                ind_tot = 1:length(vect1); ind_tot(ind_equal) = zeros(1,length(ind_equal));
                ind_diff = find(ind_tot ~= 0);
            else  % if ~isempty(ind_equal_tmp)
                ind_equal = ind_equal_tmp;
                ind_diff = (1:length(vect1));
            end  % if ~isempty(ind_equal_tmp)

        end  % if i_vect_max == 1

    % Simple case without splitting = core of the function
    else  % if (length(vect1) * length(vect2)) > mat_lim_size

        if nargin == 3

            [bid, ind_equal_tmp] = find(abs(repmat(vect1, 1, length(vect2))' - repmat(vect2,  1, length(vect1))) <= seuil_distance);
            ind_equal = vectcmp((1:max(ind_equal_tmp))',sort(ind_equal_tmp));

            ind_tot = 1:length(vect1); ind_tot(ind_equal) = zeros(1,length(ind_equal));
            ind_diff = find(ind_tot ~= 0);

        else  % if nargin == 3
            ind_equal = find(sum(repmat(vect1, 1, length(vect2))' ==...
                repmat(vect2,  1, length(vect1))));

            ind_diff = find(~(sum(repmat(vect1, 1, length(vect2))' ==...
                repmat(vect2,  1, length(vect1)))));
        end  % if nargin == 3


    end  % if (length(vect1) * length(vect2)) > mat_lim_size

end  % if (length(vect1) == 1)

return



