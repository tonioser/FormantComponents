function [pred, coef, res, moy_tot, lambda, vp_dir] = acp(Var, nb_cmp, mode, ind_var);
% [scores, basisVectors, ptsRes, meanPts, lambda] = acp(pts, nbCmp, mode, indPts);
% 
% Principal Component Analysis
% 
% Inputs
%   pts(nbObs, nbPts) : Input data points
%   nbCmp            : Number of components [default = [min(nbPts,  nbObs])]]
%   mode             : Normalisation ('norme') or not ('brut') of the scores and basis vectors [default = 'brut'] 
%   indPts           : Indices of the points considered for the calculation of the scores [default = all points]
%  
% Outpts
%   scores(nbObs, nbCmp): Scores for each basis vector, sorted in decreasing importance 
%   basisVectors(nbCmp, nbPts): Basis vectors for ALL the points
%   ptsRes(nbObs,nbPts) : Residue for ALL centred points
%   meanPts(1,nbPts)     : Mean of all points
%   lambda(nbCmp)              : Eigenvalues in decreasing importance
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 25/06/2024

%--------------------------------------------------
% Inputs and dimensions

% All observations
ind_obs = 1:size(Var,1);

% Input arguments
if ~exist('nb_cmp') nb_cmp =  min([size(Var,2), size(Var,1)]); end
if ~exist('ind_var') ind_var =  1:size(Var,2); end
if ~exist('mode') mode = 'brut'; end

% Data centring
moy_tot = mean(Var);
CVar_tot = Var - ones(size(Var,1), 1) * moy_tot;

% Keep only the relevant data points
CVar = CVar_tot(ind_obs, ind_var);
if any(isnan(CVar)) error('Some data are missing ...'); end

%--------------------------------------------------
% Eigenvalues and eigenvectors

if length(ind_var) > length(ind_obs) % Calculate the eigenvectors on the transpose
      [vp_inv, val] = eigsrt(CVar * CVar');
  lambda = diag(val);
  nb_vp = size(lambda, 1);
  vp_dir = [];
  % Recalculation of the eigenvectors for the initial matrix
  for j = 1:nb_vp
    vp_dir(:, j) = (1/sqrt(lambda(j))) * (CVar' * vp_inv(:,j));
  end
else
  [vp_dir, val] = eigsrt(CVar' * CVar);
  lambda = diag(val);
  nb_vp = size(lambda, 1);
end

%--------------------------------------------------
% Sorting, ordering, selecting

% Sort in decreasing order
lambda = lambda(length(lambda):-1:1);
vp_dir = vp_dir(:, length(lambda):-1:1);

% Keep on the nbCmp first components
ind_val = 1:nb_cmp;
lambda_min = lambda(ind_val);
vp_dir_min = vp_dir(:, ind_val);

%--------------------------------------------------
% Scores and basis vectors
% Scores
pred = CVar * vp_dir_min;
if strcmp(mode, 'norme')
  pred = pred ./ (ones(size(pred,1), 1) * std(pred));
end

% Basis vectors for all points
coef =  pred \ CVar_tot;

% ---------------------------
% Residue

CVar_tot_pred = pred * coef;
res = CVar_tot - CVar_tot_pred;

return



