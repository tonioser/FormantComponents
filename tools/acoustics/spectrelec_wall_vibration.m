function [H,YE]=spectrelec_wall_vibration(w,A,zr,l,no_vibr_paroi);

% function [H,Y] = spectrelec_wall_vibration(w,A,Zr,L,noWall);
% 
% Returns the acoutic transfer function given an area function in input,
% based on electrical analogy
% 
% Inputs
%   w(1,nbFreq)  : List of considered angular frequencies in rad.s-1
%   A(nbTubes,1) : Areas of the tubes from the glottis to the lips in cm2
%   Zr(1,nbFreq) : Radiation impedance
%   L(nbTubes,1) : Lengths of the tubes from the glottis to the lips in cm
%   noWall(1)    : Discarding from the calculation the wall vibration (1) or not (0) [default = 0]
% 
% Outputs
%   H(1,nbFreq) : Acoustic transfer function
%   Y(1,nbFreq) : Admittance
%
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 26/06/2024

%==================================================================
% Inputs & constants

c = 35900 ; % cm/s
rho = 1.098e-3 ; % g/cm^3
eta = 1.396 ;
lambda = 5.5e-5 ;
mu = 1.86e-4 ;
cp = 0.24 ;

% For wall vibration
bp = 1600 ; % dynes.s/cm
mp = 1.4 ; % g/cm²

CONST_DAT=[c ; rho ; lambda ; eta ; mu ; cp ; bp ; mp];

% Set zr if necessary
if isnan(zr(1))
    zr =  rho/(2*pi*c)*(w.^2)+j*8*rho/(3*pi.*sqrt(pi*A(end)))*w;
end  % if isnan(zr(1))

% By default take into account wall vibration
if nargin < 5
	no_vibr_paroi = 0;
end  % if nargin < 5

% All lengths equal to 1 by default
if (nargin < 4 )
   l = ones(size(A));
end;
 
c = CONST_DAT(1);
ro = CONST_DAT(2);
lambda = CONST_DAT(3);
eta = CONST_DAT(4);
mu = CONST_DAT(5);
cp = CONST_DAT(6);
bp = CONST_DAT(7);
mp = CONST_DAT(8);


%==================================================================
% Quandripoles for all tubes

%-------------------------------
% Determination of L,C,R,G,et YP

% Tube perimeters
S = 2*sqrt(A*pi) ;

% Tube characteristics
L = ro./A.*l ;
C = A.*l/ro/c/c; % (A*l)/(ro*c^2)

% Losses of Fant 1960
R_coef = sqrt(ro*mu/2*w) ;
G_coef = (eta-1)/(ro*c^2)*sqrt(lambda*w/(2*cp*ro));
R = S.*l./(A.^2) * R_coef ;
G = S.*l * G_coef ;

% Loses per wall vibration
YP_coef = 1./(bp^2+mp^2*w.^2) ;
YP = S.*l * ( (bp-j*mp*w).*YP_coef )  ;
if no_vibr_paroi; YP = 0 ; end;

% Without losses
% R = 0 ; 
% G = 0 ;
% YP = 0 ;
	
%----------------------
% Determination of Y, Z
Z = R+j*L*w ;
Y = G+j*C*w+YP ;

%---------------------------------
% Quadripole matrices
aa = (1 + (Z.*Y/2)); 
bb = - (Z + Z.^2 .* Y/4); 
cc = - Y; 
dd = aa;


%==================================================================
% Transfer function

% H=[aa bb ;
%    cc dd]

%    | A  B |
%H = |      |
%    | C  D |
%
%    | a  b |	| A  B |
%H = |      | * |      |
%    | c  d |	| C  D |

aaa = aa(1,:) ;
bbb = bb(1,:) ;
ccc = cc(1,:) ;
ddd = dd(1,:) ;

for ind = (1:length(A)-1)
	
	proda = aa(ind+1,:).*aaa + bb(ind+1,:).*ccc ;
	prodb = aa(ind+1,:).*bbb + bb(ind+1,:).*ddd ;
	prodc = cc(ind+1,:).*aaa + dd(ind+1,:).*ccc ;
	prodd = cc(ind+1,:).*bbb + dd(ind+1,:).*ddd ;
	aaa = proda ;
	bbb = prodb ;
	ccc = prodc ;
	ddd = prodd ;
end

if (zr==inf)
	H=0;
	YE= -ccc./ddd ;
else

	% Transfer function calculation
	% pl = aaa * pg + bbb * ug
	% ul = ccc * pg + ddd * ug
	% considering that pl = zr * ul
	% we get, (aaa*ddd-bbb*ccc =1)
	% ul / ug = 1 / (aaa - ccc * zr)

	H = ones(size(aaa))./( aaa - ccc.*zr) ;

	% Input conductance
	YE = -( aaa - ccc.*zr ) ./ ( bbb - ddd.*zr ) ;
end;
