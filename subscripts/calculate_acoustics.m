% Simulation of the plane wave propagation and extraction of the formants
% Takes about 4 mins for me
% 
% Author: Antoine Serrurier
% Date: 26/06/2024

%% ======================================================================
% Acoustic propagation simulation

% Overall number of articulations to treat
nbLoops = 0;
for iSpeaker = 1:length(SPEAKER)
    nbLoops = nbLoops  + size(SPEAKER(iSpeaker).cnts,1);
end  % for iSpeaker = 1:length(SPEAKER)

% Loop on the articulations to treat
hw = waitbar(0,'Simulating plane wave propagation...');
iLoop = 0;
for iSpeaker = 1:length(SPEAKER)
    for iConf = 1:size(SPEAKER(iSpeaker).cnts,1)
        
        iLoop = iLoop + 1;
        waitbar(iLoop / nbLoops,hw);

        % Simulation
        cnt = squeeze(SPEAKER(iSpeaker).cnts(iConf,:,:));
        sex = SPEAKER(iSpeaker).sex;
        [H, f, F, AF, SF, Long, VCTR, GRDVctr, midLinegrd, indGrdReg] = simulate_acoustics(...
            cnt, REG.nameArticulators, REG.indStartArticulators, REG.indEndArticulators, REG.indLandmarks, REG.nameLandmarks, sex, 0.5, constrMin);
        
        % Save
        SPEAKER(iSpeaker).Acous(iConf).H = H;
        SPEAKER(iSpeaker).Acous(iConf).f = f;
        SPEAKER(iSpeaker).Acous(iConf).F = F;
        SPEAKER(iSpeaker).Acous(iConf).AF = AF;
        SPEAKER(iSpeaker).Acous(iConf).SF = SF;
        SPEAKER(iSpeaker).Acous(iConf).Long = Long;
        SPEAKER(iSpeaker).Acous(iConf).VCTR = VCTR;
        SPEAKER(iSpeaker).Acous(iConf).GRDVctr = GRDVctr;
        SPEAKER(iSpeaker).Acous(iConf).midLinegrd = midLinegrd;
        SPEAKER(iSpeaker).Acous(iConf).indGrdReg = indGrdReg;

    end  % for iConf = 1:size(SPEAKER(iSpeaker).cnts,1)
end  % for iSpeaker = 1:length(SPEAKER)
close(hw)

%% ======================================================================
% Formants and delta-formants

for iSpeaker = 1:length(SPEAKER)
    F1 = [];
    F2 = [];
    F3 = [];
    for iConf = 1:size(SPEAKER(iSpeaker).cnts,1)
        F1 = [F1, SPEAKER(iSpeaker).Acous(iConf).F(1,2)];
        F2 = [F2, SPEAKER(iSpeaker).Acous(iConf).F(2,2)];
        F3 = [F3, SPEAKER(iSpeaker).Acous(iConf).F(3,2)];
    end  % for iConf = 1:size(SPEAKER(iSpeaker).cnts,1)
    SPEAKER(iSpeaker).F1 = F1;
    SPEAKER(iSpeaker).F2 = F2;
    SPEAKER(iSpeaker).F3 = F3;
    SPEAKER(iSpeaker).deltaF1F2 = F2-F1;
    SPEAKER(iSpeaker).deltaF2F3 = F3-F2;
end  % for iSpeaker = 1:length(SPEAKER)



