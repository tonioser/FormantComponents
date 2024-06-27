% C2 = smoothcurv(C1,P) ;
%
% Curve smoothing
% 
% Inputs
%   C1(nbPts,iXYZ): 2D or 3D input curve
%   P(1)          : filter power [default = 1]
%
% Output
%   C2(nbPts,iXYZ): 2D or 3D smoothed curve
%
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 24/06/2024

function C2 = smoothcurv(C1,P) ;

% Input arguments
if nargin == 1 ;
    P = 1 ;
end

% 2. Filter
filtre = ones(2*P+1,1) / (2*P + 1) ;

% 3. Convolution
Temp = conv2(C1,filtre) ;

% 4. Output formatting
Sup = C1(1,:) ;   % ---> Superior points
Inf = C1(end,:) ; % ---> Inferior points
for t = 1:P-1 ;
    Sup = [Sup;C1(1,:) + t*(Temp(1+2*P,:) - C1(1,:)) / P] ;
    Inf = [Temp(end-2*P,:) + t*(C1(end,:)-Temp(end-2*P,:)) / P;Inf] ;
end
C2 = [Sup;Temp(1+2*P:end-2*P,:);Inf] ;
