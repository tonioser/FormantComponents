% Plot F1 - JH relationship
% 
% Author: Antoine Serrurier
% Date: 26/06/2024


%----------------------------------------------------------------------
% Histogram of the correlation between F1 and JH for the speakers

%********************************************************
hFig = figure;
FIG2
hh = histogram(cF1JH,7);
hh.FaceAlpha = 1;
xlabelf('Pearson Coef. F1-JY');
ylabelf('Number of speakers');
set(gca, 'FontSize', 14);
%********************************************************

%-------------------------------------------------------------
% MPA morphology feature vs. F1-JH correlation

mdl = fitlm(cF1JH, MPA);
cF1JHMPACgLine = predict(mdl, cF1JH);

%***********************************************************
hFig = figure;
FIG2
plot(cF1JH, MPA, '.', 'MarkerSize', 18)
plot(cF1JH, cF1JHMPACgLine, 'LineWidth', 4);
xlabelf('Correlation coef. F1-JY');
ylabelf('MPA');
set(gca, 'FontSize', 14);
%***********************************************************


