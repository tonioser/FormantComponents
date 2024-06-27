% [V, D] = eigsrt(X)
% 
% Returns the eigenvectors and eigenvalues of a matrix X
% 

function [V, D] = eigsrt(X)
[V, D] = eig(X);
[vp_ord ind_ord] = sort(diag(D));
D = diag(vp_ord);
V = V(:,ind_ord);
return
