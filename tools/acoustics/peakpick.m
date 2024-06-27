% [maxi, ind_maxi, mini, ind_mini] = peakpick(signal);
% 
% Returns the local maximums and minimums of a signal
% 
% Inputs
%   signal(1:nb_ech) : input signal
% 
% Outputs
%   maxi, ind_maxi   : values and indices of the local maximums
%   mini, ind_mini   : values and indices of the local minimums
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 26/06/2024

function [maxi, ind_maxi, mini, ind_mini] = peakpick(signal);

% Initialisation
mini = []; ind_mini = []; maxi = []; ind_maxi = [];

nb_ech = length(signal);
for ind = 2:nb_ech-1
  if (signal(ind-1) < signal(ind)) & (signal(ind+1) < signal(ind))
    ind_maxi = [ind_maxi, ind]; maxi = [maxi, signal(ind)];
  end
  if (signal(ind-1) > signal(ind)) & (signal(ind+1) > signal(ind))
    ind_mini = [ind_mini, ind]; mini = [mini, signal(ind)];
  end
end

%  if no local extremum, returns the signal extremty
if isempty(maxi) 
	[maxi, ind_maxi] = max(signal);
end
if isempty(mini) 
	[mini, ind_mini] = min(signal);
end
return
