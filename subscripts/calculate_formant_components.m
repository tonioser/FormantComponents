% Calculate the formant and delta-formant components
% 
% Author: Antoine Serrurier
% Date: 26/06/2024

nameFmts = {'F1', 'F2', 'F3', 'deltaF1F2', 'deltaF2F3'};
nbFmts = length(nameFmts);

%------------------------------------------------------------
% Multiple linear regression of the vowel cnts on the formants

for iFmt = 1:nbFmts
    namePred = nameFmts{iFmt};
    for iSpeaker = 1:length(SPEAKER)
        
        % Predictor = formant value
        eval(['ZP = SPEAKER(iSpeaker).', namePred, ''';']);
        % Z-scored
        ZP = (ZP - mean(ZP)) / std(ZP);
        
        % Variables
        pts = SPEAKER(iSpeaker).cnts;
        
        % Multiple linear regressions
        [basisVectorRaw, cntMean, res] = regress_Data_Scores_2_BasisVector_adapted(pts, ZP);
        basisVector = NaN([1, size(basisVectorRaw)]);
        basisVector(1,:,:) = basisVectorRaw;
        
        % Save
        SPEAKER(iSpeaker).RL(iFmt).ZP = ZP;
        SPEAKER(iSpeaker).RL(iFmt).pts = pts;
        SPEAKER(iSpeaker).RL(iFmt).basisVector = basisVector;
        SPEAKER(iSpeaker).RL(iFmt).res = res;
        SPEAKER(iSpeaker).RL(iFmt).cntMean = cntMean;

    end  % for iSpeaker = 1:length(SPEAKER)
end  % for iFormant = 1:2

%------------------------------------------------------------
% Gather all basis vectors and cntMean

basisVectors = NaN(nbFmts, nbSpeakers, size(SPEAKER(1).RL(1).basisVector,2), 2);
cntMeans = NaN(nbFmts, nbSpeakers, size(SPEAKER(1).RL(1).cntMean,1), 2);
for iFmt = 1:length(nameFmts)
    for iSpeaker = 1:length(SPEAKER)
        basisVectors(iFmt,iSpeaker,:,:) = SPEAKER(iSpeaker).RL(iFmt).basisVector;
        cntMeans(iFmt,iSpeaker,:,:) = SPEAKER(iSpeaker).RL(iFmt).cntMean;
    end  % for iSpeaker = 1:length(SPEAKER)
end  % for iFormant = 1:2

%------------------------------------------------------------
% Extremum values of the preds

ZPsExtr = NaN(nbFmts, nbSpeakers, 2);
for iFmt = 1:length(nameFmts)
    for iSpeaker = 1:length(SPEAKER)
        ZPsExtr(iFmt,iSpeaker,1) = min(SPEAKER(iSpeaker).RL(iFmt).ZP);
        ZPsExtr(iFmt,iSpeaker,2) = max(SPEAKER(iSpeaker).RL(iFmt).ZP);
    end  % for iSpeaker = 1:length(SPEAKER)
end  % for iFormant = 1:2

