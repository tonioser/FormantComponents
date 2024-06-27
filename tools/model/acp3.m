function [pred, coef, res, moy, lambda, vp_dir] = acp3(Pts, nb_cmp, norm_mode, ind_var);
% [scores, basisVectors, ptsRes, meanPts, lambda] = acp3(pts, nbCmp, mode);
% 
% Principal Component Analysis (including 3D data)
% 
% Inputs
%   pts(nbObs, nbPts, iXYZ) : Input data points
%   nbCmp                   : Number of components [default = [min(nbPts,  nbObs])]]
%   mode                    : Normalisation ('norme') or not ('brut') of the scores and basis vectors [default = 'brut'] 
%  
% Outpts
%   scores(nbObs, nbCmp)             : Scores for each basis vector, sorted in decreasing importance 
%   basisVectors(nbCmp, nbPts, iXYZ) : Basis vectors for ALL the points
%   ptsRes(nbObs,nbPts, iXYZ)        : Residue for ALL centred points
%   meanPts(nbPts, iXYZ)             : Mean of all points
%   lambda(nbCmp)                    : Eigenvalues in decreasing importance
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 25/06/2024
%

if ~exist('norm_mode') norm_mode = 'brut'; end

% Reshaping 3D -> 1D
[nb_obs, nb_pts, ndim] = size(Pts);
Var = reshape(Pts, nb_obs, nb_pts*ndim);

% 1D PCA
[pred, coef_, res_, moy_, lambda, vp_dir] = acp(Var, nb_cmp, norm_mode);

% Reshaping 1D -> 3D
coef = reshape(coef_, nb_cmp, nb_pts, ndim);
res = reshape(res_, nb_obs, nb_pts, ndim);
moy = squeeze(reshape(moy_, 1, nb_pts, ndim));

return

