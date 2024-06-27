% Plot the correlations between inter-speaker components and morphology features
% 
% Author: Antoine Serrurier
% Date: 26/06/2024


%*************************************************************
% MPA vs. 1st inter-speaker comp. of F1

mdl = fitlm(PCA2ndLevel(1).scores(:,1),MPA);
MPACgline = predict(mdl, PCA2ndLevel(1).scores(:,1));
c = corrcoef([PCA2ndLevel(1).scores(:,1),MPA]);
c = round(c(1,2)*100) / 100;

cols = getColPlot(3);
hFig = figure;
FIG2
plot(PCA2ndLevel(1).scores(indM,1), MPA(indM), '.', 'MarkerSize', 30)
plot(PCA2ndLevel(1).scores(indF,1), MPA(indF), '.', 'MarkerSize', 30)
plot(PCA2ndLevel(1).scores(:,1), MPACgline, 'LineWidth', 4, 'Color', cols(3,:))
XL = xlim;
YL = ylim;
xlim([-2.8, XL(2)])
text(-1.5, -5.5, '1^{st} inter-speaker comp.  of F1', 'FontSize', 14)
text(-2.6, -5.3, '<- Anterior        MPA        Posterior ->', 'FontSize', 14, 'Rotation', 90)
text(-2, 5, ['R = ', num2str(c)], 'FontSize', 18, 'FontWeight', 'b', 'BackgroundColor', [0.9 0.9 0.9])
set(gca, 'FontSize', 14)
%*************************************************************



%*************************************************************
% MPC vs. 1st inter-speaker comp. of F1

mdl = fitlm(PCA2ndLevel(1).scores(:,1),MPC);
MPCgline = predict(mdl, PCA2ndLevel(1).scores(:,1));
c = corrcoef([PCA2ndLevel(1).scores(:,1),MPC]);
c = round(c(1,2)*100) / 100;

cols = getColPlot(3);
hFig = figure;
FIG2
plot(PCA2ndLevel(1).scores(indM,1), MPC(indM), '.', 'MarkerSize', 30)
plot(PCA2ndLevel(1).scores(indF,1), MPC(indF), '.', 'MarkerSize', 30)
plot(PCA2ndLevel(1).scores(:,1), MPCgline, 'LineWidth', 4, 'Color', cols(3,:))
XL = xlim;
YL = ylim;
xlim([-2.8, XL(2)])
text(-1.5, 2.12, '1^{st} inter-speaker comp.  of F1', 'FontSize', 14)
text(-2.6, 2.2, '<- Domed        MPC (cm)        Flat ->', 'FontSize', 14, 'Rotation', 90)
text(-2, 4, ['R = ', num2str(c)], 'FontSize', 18, 'FontWeight', 'b', 'BackgroundColor', [0.9 0.9 0.9])
set(gca, 'FontSize', 14)
%*************************************************************


%*************************************************************
% MY vs. 1st inter-speaker comp. of F2

mdl = fitlm(PCA2ndLevel(2).scores(:,1),MY);
MYgline = predict(mdl, PCA2ndLevel(2).scores(:,1));
c = corrcoef([PCA2ndLevel(2).scores(:,1),MY]);
c = round(c(1,2)*100) / 100;

hFig = figure;
FIG2
plot(PCA2ndLevel(2).scores(indM,1), 10 - MY(indM), '.', 'MarkerSize', 30)
plot(PCA2ndLevel(2).scores(indF,1), 10 - MY(indF), '.', 'MarkerSize', 30)
plot(PCA2ndLevel(2).scores(:,1), 10-MYgline, 'LineWidth', 4)
set(gca, 'YDir', 'Reverse')
XL = xlim;
xlim([-2.2, XL(2)])
text(-1, 7.8, '1^{st} inter-speaker comp. of F2', 'FontSize', 14, 'BackgroundColor', 'w')
text(-2.05, 7.9, '<- Longer       MY (cm)       Shorter ->', 'FontSize', 14, 'Rotation', 90)
text(-1.5, 5.3, ['R = ', num2str(c)], 'FontSize', 18, 'FontWeight', 'b', 'BackgroundColor', [0.9 0.9 0.9])
set(gca, 'FontSize', 14)
%*************************************************************


%*************************************************************
% MY vs. 1st inter-speaker comp. of deltaF1F2

mdl = fitlm(PCA2ndLevel(4).scores(:,1),MY);
MYgline = predict(mdl, PCA2ndLevel(4).scores(:,1));
c = corrcoef([PCA2ndLevel(4).scores(:,1),MY]);
c = round(c(1,2)*100) / 100;

hFig = figure;
FIG2
plot(PCA2ndLevel(4).scores(indM,1), 10 - MY(indM), '.', 'MarkerSize', 30)
plot(PCA2ndLevel(4).scores(indF,1), 10 - MY(indF), '.', 'MarkerSize', 30)
plot(PCA2ndLevel(4).scores(:,1), 10-MYgline, 'LineWidth', 4)
set(gca, 'YDir', 'Reverse')
XL = xlim;
xlim([-2.5, XL(2)])
text(-1.5, 7.8, '1^{st} inter-speaker comp. of \DeltaF1F2', 'FontSize', 14, 'BackgroundColor', 'w')
text(-2.35, 7.9, '<- Longer       MY (cm)       Shorter ->', 'FontSize', 14, 'Rotation', 90)
text(-1.8, 5.3, ['R = ', num2str(c)], 'FontSize', 18, 'FontWeight', 'b', 'BackgroundColor', [0.9 0.9 0.9])
set(gca, 'FontSize', 14)
%*************************************************************


%*************************************************************
% MY vs. 1st inter-speaker comp. of deltaF2F3

mdl = fitlm(PCA2ndLevel(5).scores(:,1),MY);
MYgline = predict(mdl, PCA2ndLevel(5).scores(:,1));
c = corrcoef([PCA2ndLevel(5).scores(:,1),MY]);
c = round(c(1,2)*100) / 100;

hFig = figure;
FIG2
plot(PCA2ndLevel(5).scores(indM,1), 10 - MY(indM), '.', 'MarkerSize', 30)
plot(PCA2ndLevel(5).scores(indF,1), 10 - MY(indF), '.', 'MarkerSize', 30)
plot(PCA2ndLevel(5).scores(:,1), 10-MYgline, 'LineWidth', 4)
set(gca, 'YDir', 'Reverse')
XL = xlim;
xlim([-2.2, XL(2)])
text(-1.2, 7.8, '1^{st} inter-speaker comp. of \DeltaF2F3', 'FontSize', 14, 'BackgroundColor', 'w')
text(-2.1, 7.9, '<- Longer       MY (cm)       Shorter ->', 'FontSize', 14, 'Rotation', 90)
text(-1.5, 5.3, ['R = ', num2str(c)], 'FontSize', 18, 'FontWeight', 'b', 'BackgroundColor', [0.9 0.9 0.9])
set(gca, 'FontSize', 14)
%*************************************************************


