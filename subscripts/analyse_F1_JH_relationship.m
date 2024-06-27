% Analyse F1 - JH relationship
% 
% Author: Antoine Serrurier
% Date: 26/06/2024


%------------------------------------------------------------
% Calculate JH = highest points of the lower incisors

iJLT = REG.indLandmarks(find(strcmp(REG.nameLandmarks, 'JLT')));

for iSpeaker = 1:length(SPEAKER)
    JH = SPEAKER(iSpeaker).cnts(:,iJLT,2);
    ZJH = (JH - mean(JH)) / std(JH);
    SPEAKER(iSpeaker).JH = JH;
    SPEAKER(iSpeaker).ZJH = ZJH;
end  % for iSpeaker = 1:length(SPEAKER)

%-------------------------------------------------------------
% Gather F1 and JH

% Overall
F1s = [];
JHs = [];
for iSpeaker = 1:length(SPEAKER)
    F1s = [F1s, SPEAKER(iSpeaker).F1];
    JHs = [JHs, SPEAKER(iSpeaker).JH'];
end  % for iSpeaker = 1:length(SPEAKER)

%-------------------------------------------------------------
% Correlation F1-JH for all speakers

% Loops on the speakers
cF1JH = [];
for iSpeaker = 1:length(SPEAKER)
    c = corrcoef([SPEAKER(iSpeaker).F1', SPEAKER(iSpeaker).JH]);
    cF1JH = [cF1JH, c(1,2)];
end  % for iSpeaker = 1:length(SPEAKER)
cF1JH = abs(cF1JH)';

mean(cF1JH);
disp(['Mean correlation F1-JH across speakers = ', num2str(mean(cF1JH))])


