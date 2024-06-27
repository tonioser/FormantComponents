% Plot second-level component nomograms
% 
% Author: Antoine Serrurier
% Date: 26/06/2024

%--------------------------------------------------------------------------
% Overall mean articulation of the data
cntsMean = [];
for iSpeaker = 1:nbSpeakers
    cntsMean = cat(3,cntsMean,SPEAKER(iSpeaker).RL(1).cntMean);
end  % for iSpeaker = 1:nbSpeakers
cntMean = squeeze(mean(cntsMean,3));


%--------------------------------------------------------------------------
% Figures = loop on the formants

iScores2ndLevel = 1
hFigs = NaN(3,2);
XL = [];
YL = [];
for iFmt = 1:nbFmts
    for iExtr = 1:2
        switch iExtr
            case 1  % switch iExtr
                scoresVal = min(PCA2ndLevel(iFmt).scores(:,iScores2ndLevel));
            case 2  % switch iExtr
                scoresVal = max(PCA2ndLevel(iFmt).scores(:,iScores2ndLevel));
        end  % switch iExtr

        % Predict the first level basis vectors for the two extreme cases
        ptsEst = squeeze(predict_Scores_BasisVectors_2_Data(scoresVal, PCA2ndLevel(iFmt).basisVector(iScores2ndLevel,:,:), PCA2ndLevel(iFmt).moy));
        basisVectorEst = ptsEst;

        % Nomograms
        [cntNomos, scoresNomos] = nomos2D(cntMean, ZPExtr(iFmt,:)', basisVectorEst);

        %*************************************************************
        iScoreNomo = 1;
        hFigs(iFmt, iExtr) = figure;
        FIG
        for iNomo = 1:size(cntNomos,2)
            if scoresNomos(iScoreNomo, iNomo) <= 0
                c = colR;
            else  % if predNomos(iPred, iNomo) <= 0
                c = colG;
            end  % if predNomos(iPred, iNomo) <= 0
            cnt = squeeze(cntNomos(iScoreNomo,iNomo,:,:));
            cnt(REG.indC5,:) = [];
            plotf(cnt, 'Color', c, 'LineWidth', 2);
            plotf(cnt(1:20:end,:), 'k.', 'MarkerSize', 5);
        end  % for iNomo = 1:size(cntNomos,2)
        cntMeanPlot = cntMean;
        cntMeanPlot(REG.indC5,:) = [];
        plotf(cntMeanPlot, 'k', 'LineWidth', 3);
        XL = [XL; xlim];
        YL = [YL; ylim];
        set(gca, 'FontSize', 14)
        titlef([nameFmts{iFmt}, ' comp. - 2nd-level comp.', num2str(iScores2ndLevel), ' Extr ', num2str(iExtr)]);
        %*************************************************************
    end  % for iExtr = 1:2
end  % for iFmt = 1 = [1,2,3]
XLm = [min(XL(:,1)), max(XL(:,2))];
YLm = [min(YL(:,1)), max(YL(:,2))];
for iFmt = [1,2,3,4,5]
    for iExtr = 1:2
        %*************************************************************
        figure(hFigs(iFmt, iExtr))
        axis([XLm, YLm])
        text(8, 3.5, 'X (cm)', 'FontSize', 14, 'BackgroundColor', 'w')
        text(2, 8, 'Y (cm)', 'FontSize', 14, 'BackgroundColor', 'w', 'Rotation', 90)
        xlabel('')
        ylabel('')
        drawnow
        %*************************************************************
    end  % for iExtr = 1:2
end  % for iFmt = 1 = [1,2,3]

