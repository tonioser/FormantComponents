% Analysis and inter-speaker analysis of the formant components
% 
% Replicate the code described in the following article:
% Antoine Serrurier and Christiane Neuschaefer-Rube (2024, in review)
% Formant-based articulatory strategies: characterisation and inter-speaker variability analysis
% Journal of Phonetics
%
% The subscripts need to be called one after the other: the place is not flexible. 
% The subscripts calling figures are optional, but their place is not flexible either. 
% 
% The code does as follows:
%   - Set the path
%   - Load constants and data
%   - Calculate the morphology features
%   - Simulate plane acoustic wave propagation and extract the formants
%   - Calculate the formant and delta-formant components
%   - Analyse the F1-JH relationships
%   - Calculate the second-level inter-speaker components
%   - Analyse the relationship between the inter-speaker components and the morphology features 
% 
% Any use and this code should cite the following paper:
% Antoine Serrurier and Christiane Neuschaefer-Rube (2024, in review)
% Formant-based articulatory strategies: characterisation and inter-speaker variability analysis
% Journal of Phonetics
% 
% Author: Antoine Serrurier
% Date: 26/06/2024

%--------------------------------------------------------------------------
% Set path
addpath(genpath('./subscripts/'))
addpath(genpath('./tools/'))

%--------------------------------------------------------------------------
% Load constants and data
set_constants_load_data

%--------------------------------------------------------------------------
% Calculate the morphology features
calculate_morphology_features

%--------------------------------------------------------------------------
% Simulate plane acoustic wave propagation and extract the formants
calculate_acoustics

%--------------------------------------------------------------------------
% Calculate the formant and delta-formant components
calculate_formant_components
% Plot some cases (optional)
plot_formant_components_samples
% Plot averaged formant components (optional)
plot_average_formant_components

%--------------------------------------------------------------------------
% Analyse the F1-JH relationships
analyse_F1_JH_relationship
% Plot the F1-JH relationships (optional)
plot_F1_JH_relationship

%--------------------------------------------------------------------------
% Calculate inter-speaker components
calculate_secondlevel_components
% Plot inter-speaker components (optional)
plot_secondlevel_components

%--------------------------------------------------------------------------
% Analyse the relationship between the inter-speaker components and the morphology 
analyse_secondlevelcomp_morphology_relationship
% Plot the correlations between inter-speaker components and morphology features (optional)
plot_secondlevelcomp_morphology_relationship



