function [MPC, MPCC] = gPCA_getMPC_adapted(cnts, indPalCurv)
% 
% Calculate the MPC control parameter for a set of articulation contours.
%
% Inputs
%     cnts(nbSubjects,nbPts,nbDim) : Articulation contours
%                                    Typically of size 41 x 1692 x 2
%     indPal(nbPtsPal)             : Indices of the vocal tract points corresponding to the curved part of the hard palate for an articulation contour
%                                    Typically of length 69
% 
% Outputs
%     MPC(nbSubjects)                 : Column vector of the control parameter MPC
%     MPCC(nbSubjects)                : Column vector of the centred control parameter MPC
% 
% Author : Antoine Serrurier
% Date: 21/01/2024

% Sizes
nbObs = size(cnts,1);

% Calculate the least-square circle on the curved points of the palate
% MPC defined as its radius
MPC = [];
for iSpeaker = 1:nbObs
    cnt = squeeze(cnts(iSpeaker, indPalCurv, :));
    circle = circle_least_square(cnt);
    MPC = [MPC; circle.radius];
end  % for iSpeaker = 1:nbObs

% Centred
MPC;
MPCC = MPC - mean(MPC);

end