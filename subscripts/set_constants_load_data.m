% Define the constants and load the data
% 
% Author: Antoine Serrurier
% Date: 26/06/2024

%% ======================================================================
% Set constants

constrMin = 0.1; % Imposed minimal constriction of the rea function (in cm2) to calculate the acoustic transfer function

% Green and red colors forthe plots
colG = [0.4660    0.6740    0.1880];
colR = [0.8500    0.3250    0.0980];

%% ======================================================================
% Load data

dirData = './data';

%-------------------------------------------------
% Speakers

% Contours
load([dirData, '/SPEAKERS.mat']);

% Number of speakers
nbSpeakers = length(SPEAKER);

% Indices of the females and males speakers
indF = find([SPEAKER.sex] == 'f');
indM = find([SPEAKER.sex] == 'm');

%-------------------------------------------------
% Regions and landmarks

load([dirData, '/REGIONS.mat']);

% REG.UT = index of the lowest point of the upper incisors
% 
% REG.nameArticulators     = names of the articulators in the contours
% REG.indStartArticulators = indices of the start of each articulator in the contours
% REG.indEndArticulators   = indices of the end of each articulator in the contours
% REG.indLandmarks         = Indices of the anatomical landmarks in the contours
% REG.nameLandmarks        = Names of the anatomical landmarks
% 
% REG.indOrgsVT     = Indices in the contours corresponding to the vocal tract points 
% REG.indPalVT      = Indices in the contours corresponding to the vocal tract points of the hard palate 
% REG.indPalVT_curv = Indices in REG.indPalVT corresponding to the curved part of the palate
% REG.indPhaVT      = Indices in the contours corresponding to the vocal tract points of the pharyngeal wall 
% REG.indC5         = Indices in the contours corresponding to the C5 verterbrae 

