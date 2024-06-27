% [nb_zones, ind_start, ind_end] = nb_contig_zones(X);
% 
% Returns the number of zones of a vector split by NaN
% 
% Inputs
%   X(1,nbPts): Input vector
% 
% Outputs
%   nb_zones(1)            : Number of zones
%   ind_start(1, nb_zones) : Indices in X of the beginning of the zones
%   ind_end(1, nb_zones)   : Indices in X of the end of the zones
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 24/06/2024


function [nb_zones, ind_start, ind_end] = nb_contig_zones(Xcnt);

if ~isempty(Xcnt)
    [NbCont, NbPts] = size(Xcnt);
    for NoCont = 1:NbCont
        ind_tmp = find(~isnan(Xcnt(NoCont,:)));
        if ~isempty(ind_tmp)
            tmp = [Xcnt(NoCont,ind_tmp(1):NbPts), NaN]; % Just add a NaN at the end
            tmp_AnA = ~isnan(tmp);
            tmp_absdiff = abs(diff(tmp_AnA));
            nb_zones(NoCont,1) = ceil((sum(tmp_absdiff)+1)/2);
            tmp_diff = diff(tmp_AnA);
            % First index of start of first zone
            ind_start(NoCont, 1) = ind_tmp(1);
            ind_fin_tmp = find(tmp_diff == -1);
            ind_end(NoCont, 1:length(ind_fin_tmp)) = ind_fin_tmp + ind_tmp(1) -1;
            ind_start(NoCont, 1) = ind_tmp(1);
            ind_deb_tmp = find(tmp_diff == 1);
            ind_start(NoCont, 2:length(ind_deb_tmp)+1) = ind_deb_tmp + ind_tmp(1);
        else
            nb_zones(NoCont,1) = 0; ind_start = 0; ind_end = 0;
        end
    end

    ind_start = repl_val(ind_start, 0, NaN);
    ind_end = repl_val(ind_end, 0, NaN);

else
    nb_zones = 0; ind_start = 0; ind_end = 0;
end
return
