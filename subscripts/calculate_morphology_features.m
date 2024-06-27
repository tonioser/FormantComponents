% Calculate the morphology features
% 
% The code in this section is based on the methods detailed in:
% A. Serrurier, C. Neuschaefer-Rube
% Morphological and acoustic modeling of the vocal tract
% The Journal of the Acoustical Society of America, 153 (3) (2023) 1867–1886
% doi:10.1121/10.0017356
% 
% Author: Antoine Serrurier
% Date: 26/06/2024

% Gather the average-articulations of all speakers
cntsMorph = [];
for iSpeaker = 1:length(SPEAKER)
    cntsMorph = cat(3, cntsMorph, SPEAKER(iSpeaker).morph);
end  % for iSpeaker = 1:length(SPEAKER)
cntsMorph = permute(cntsMorph, [3,1,2]);

% Average average-articulation
meanMorph = squeeze(mean(cntsMorph));

% Replicate the average average-articulation to be of same size than cntsMorph
meanMorphRep = permute(repmat(meanMorph, [1,1,nbSpeakers]),[3,1,2]);

% MX morphological feature
[MX, MXC] = gPCA_getMX(cntsMorph, REG.UT(2), REG.indPhaVT);
% Residue
[~, ~, resMX] = regress_Data_Scores_2_BasisVector_adapted(cntsMorph, MXC);

% MY morphological feature
iGF = REG.indLandmarks(find(strcmp(REG.nameLandmarks, 'GF')));
iGB = REG.indLandmarks(find(strcmp(REG.nameLandmarks, 'GB')));
[MY, MYC] = gPCA_getMY(resMX + meanMorphRep, iGF, iGB);
% Residue
[~, ~, resMXMY] = regress_Data_Scores_2_BasisVector_adapted(resMX, MYC);

% MA morphological feature
iPhL = REG.indLandmarks(find(strcmp(REG.nameLandmarks, 'PhL')));
iPhU = REG.indLandmarks(find(strcmp(REG.nameLandmarks, 'PhU')));
[MA, MAC] = gPCA_getMA(resMXMY + meanMorphRep, iPhL, iPhU);
% Residue
[~, ~, resMXMYMA] = regress_Data_Scores_2_BasisVector_adapted(resMXMY, MAC);

% MPA morphological feature
[MPA, MPAC] = gPCA_getMPA(resMXMYMA + meanMorphRep, REG.indPalVT);
[~, ~, resMXMYMAMPA] = regress_Data_Scores_2_BasisVector_adapted(resMXMYMA, MPAC);

% MPC morphological feature
[MPC, MPCC] = gPCA_getMPC_adapted(resMXMYMAMPA + meanMorphRep, REG.indPalVT(REG.indPalVT_curv));

