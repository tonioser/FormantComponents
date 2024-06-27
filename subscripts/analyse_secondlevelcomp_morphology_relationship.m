% Analyse inter-speaker components - morphology relationships
% 
% Author: Antoine Serrurier
% Date: 26/06/2024

% Loop on the formants to calculate the correlation coefficient
Cmat = [];
for iFmt = 1:nbFmts

    % Correlation coefficients
    C = corrcoef([PCA2ndLevel(iFmt).scores, MX, MY, MA, MPA, MPC]);

    % Cleaning and formatting
    C = round(abs(C)*100) / 100;
    C(size(PCA2ndLevel(iFmt).scores,2)+1:end,:) = [];
    C(:,1:size(PCA2ndLevel(iFmt).scores,2)) = [];

    % Save
    Cmat = [Cmat; C];
end  % for iFmt = [1,2,3]

% Display
tg = array2table(Cmat);
tg.Properties.VariableNames = ["MX", "MY", "MA", "MPA", "MPC"];
tg.Properties.RowNames = [...
    "1st inter-speaker component of F1 comp.",...
    "1st inter-speaker component of F2 comp.",...
    "1st inter-speaker component of F3 comp.",...
    "1st inter-speaker component of deltaF1F2 comp.",...
    "1st inter-speaker component of deltaF2F3 comp."];
disp(tg)


