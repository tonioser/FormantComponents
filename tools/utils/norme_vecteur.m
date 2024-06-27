% function A = norme_vecteur(A);
% 
% Normalisation of the rows of a marix A
% 


function A = norme_vecteur(A)
warning off

% Norms
NA = Norm2(A) ;

% Loop
for t = 1:size(A,2) ;
    A(:,t) = A(:,t) ./ NA ;
end
