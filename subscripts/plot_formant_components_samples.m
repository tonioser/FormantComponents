% Plot formant components nomograms for particular speakers and formants
% 
% Author: Antoine Serrurier
% Date: 26/06/2024

%---------------------------------------------------------------
% Selection of the components to plot

% Formant
iFmt = 2;

% List of speakers
indSelSubs = [3, 15, 24, 32];

%---------------------------------------------------------------
hFigs = [];
XL = [];
YL = [];
for iFig = 1:length(indSelSubs)
    iSpeaker = indSelSubs(iFig);
    
    [cntNomos, predNomos] = nomos2D(SPEAKER(iSpeaker).RL(iFmt).cntMean,...
        SPEAKER(iSpeaker).RL(iFmt).ZP, SPEAKER(iSpeaker).RL(iFmt).basisVector);
    iPred = 1;
    
    %*************************************************************
    hFigs(iFig) = figure;
    FIG
    for iNomo = 1:size(cntNomos,2)
        if predNomos(iPred, iNomo) <= 0
            c = colR;
        else  % if predNomos(iPred, iNomo) <= 0
            c = colG;
        end  % if predNomos(iPred, iNomo) <= 0
        cnt = squeeze(cntNomos(iPred,iNomo,:,:));
        cnt(REG.indC5,:) = [];
        plotf(cnt, 'Color', c, 'LineWidth', 2);
        % plotf(cnt(1:20:end,:), '.', 'Color', c);
        plotf(cnt(1:20:end,:), 'k.', 'MarkerSize', 5);
    end  % for iNomo = 1:size(cntNomos,2)
    cntMeanPlot = SPEAKER(iSpeaker).RL(iFmt).cntMean;
    cntMeanPlot(REG.indC5,:) = [];
    plotf(cntMeanPlot, 'k', 'LineWidth', 3);
    titlef(['Speaker ', num2str(iSpeaker), ' - ', nameFmts{iFmt}]);
    XL = [XL; xlim];
    YL = [YL; ylim];
    set(gca, 'FontSize', 14)
    %*************************************************************
end  % for iFig = 1:length(indSelSubs)
for iFig = 1:length(hFigs)
    %*************************************************************
    figure(hFigs(iFig))
    axis([min(XL(:,1)), max(XL(:,2)), min(YL(:,1)), max(YL(:,2))])
    text(7.5, 2, 'X (cm)', 'FontSize', 16, 'BackgroundColor', 'w');
    text(1, 6, 'Y (cm)', 'FontSize', 16, 'BackgroundColor', 'w', 'Rotation', 90);
    xlabel('');
    ylabel('');
    %*************************************************************
end  % for iSex = 1:2

