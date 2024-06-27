% Calculate second-level components
% 
% Author: Antoine Serrurier
% Date: 26/06/2024

%------------------------------------------------------------
% Gather all basis vectors

basisVectors = NaN(nbFmts, nbSpeakers, size(SPEAKER(1).RL(1).basisVector,2), 2);
for iFmt = 1:length(nameFmts)
    for iSpeaker = 1:length(SPEAKER)
        basisVectors(iFmt,iSpeaker,:,:) = SPEAKER(iSpeaker).RL(iFmt).basisVector;
    end  % for iSpeaker = 1:length(SPEAKER)
end  % for iFormant = 1:2

%------------------------------------------------------------
% PCA of the basis vectors of the first level = 2nd level PCA

clear PCA2ndLevel
for iFmt = 1:nbFmts
    
    % Basis vectors of the 1st level for the VT points only
    pts = squeeze(basisVectors(iFmt,:,REG.indOrgsVT,:));
    indAnA = find(~isnan(pts(1,:,1)));
    ptsAnA = pts(:,indAnA,:);
    
    % PCA
    nbPred2ndLevel = 1;
    [scores2nd, basisVector2ndVT, ~, meanPts2nd] = acp3(ptsAnA, nbPred2ndLevel, 'norme');

    % Percentage of variance explanation
    ptsAnAEst = predict_Scores_BasisVectors_2_Data(scores2nd, basisVector2ndVT, meanPts2nd);
    varex_tot = variance_rms(ptsAnA, ptsAnAEst);
    disp(['Percentage of variance explanation ', nameFmts{iFmt}, ' = ', num2str(round(varex_tot*100)), '%'])
    
    % Linear regression to predict together the non-VT points together
    pts = squeeze(basisVectors(iFmt,:,:,:));
    [basisVector2ndRaw, meanPts2nd, resPts2nd] = regress_Data_Scores_2_BasisVector_adapted(pts, scores2nd);
    basisVector2nd = NaN([nbPred2ndLevel, size(basisVector2ndRaw)]);
    basisVector2nd(1,:,:) = basisVector2ndRaw;
    
    % Store
    PCA2ndLevel(iFmt).pts = pts;
    PCA2ndLevel(iFmt).scores = scores2nd;
    PCA2ndLevel(iFmt).basisVector = basisVector2nd;
    PCA2ndLevel(iFmt).res = resPts2nd;
    PCA2ndLevel(iFmt).moy = meanPts2nd;
    PCA2ndLevel(iFmt).varex_tot = varex_tot;

end  % for iFmt = 1:nbFmts

%------------------------------------------------------------
% Extremum values of the preds

ZPsExtr = NaN(nbFmts, nbSpeakers, 2);
for iFmt = 1:length(nameFmts)
    for iSpeaker = 1:length(SPEAKER)
        ZPsExtr(iFmt,iSpeaker,1) = min(SPEAKER(iSpeaker).RL(iFmt).ZP);
        ZPsExtr(iFmt,iSpeaker,2) = max(SPEAKER(iSpeaker).RL(iFmt).ZP);
    end  % for iSpeaker = 1:length(SPEAKER)
end  % for iFormant = 1:2
ZPExtr = squeeze(mean(ZPsExtr,2));

