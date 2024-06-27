% Plot averaged formant and delta-formant components nomograms
% 
% Author: Antoine Serrurier
% Date: 26/06/2024

hFigs = [];
XL = [];
YL = [];
iFig = 0;
for iFmt = 1:nbFmts
    iFig = iFig + 1;

    % Average basis vectors
    basisVectorMeanSex = mean(squeeze(basisVectors(iFmt,:,:,:)));
    % Average mean
    cntMeanSex = squeeze(mean(cntMeans(iFmt, :, :, :)));
    % Average extreme ZPs
    ZPsMeanSex = mean(squeeze(ZPsExtr(iFmt,:,:)));

    % Nomograms
    [cntNomos, predNomos] = nomos2D(cntMeanSex, ZPsMeanSex', basisVectorMeanSex);
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
    cntMeanPlot = cntMeanSex;
    cntMeanPlot(REG.indC5,:) = [];
    plotf(cntMeanPlot, 'k', 'LineWidth', 3);
    XL = [XL; xlim];
    YL = [YL; ylim];
    set(gca, 'FontSize', 14)
    titlef([nameFmts{iFmt}]);
    %*************************************************************

end  % for iFmt = 1:nbFmts
for iFig = 1:length(hFigs)
    %*************************************************************
    figure(hFigs(iFig))
    axis([min(XL(:,1)), max(XL(:,2)), min(YL(:,1)), max(YL(:,2))])
    text(7.5, 3.8, 'X (cm)', 'FontSize', 14, 'BackgroundColor', 'w');
    text(2, 7.5, 'Y (cm)', 'FontSize', 14, 'BackgroundColor', 'w', 'Rotation', 90);
    xlabel('');
    ylabel('');
    %*************************************************************
end  % for iSex = 1:2

