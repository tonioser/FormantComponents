function pp = csape_(x,y,conds,valconds)
%CSAPE	Cubic spline interpolation with various end-conditions.
%
%        pp = csape(x,y[,conds[,valconds]])
%
% returns the cubic spline interpolant (in pp-form) to the given data  (x,y)
% using the specified end-conditions  conds(i)  with  values  valconds(i) ,
% with i=1 (i=2) referring to the left (right) endpoint. 
%   conds(i)=j  means that the j-th derivative is being specified to be
%   valconds(i) , j=1,2. 
%   conds(1)=0=conds(2)  means periodic end conditions.
%  If conds(i) is not specified or is different from 0, 1 or 2, then the 
%  default value for  conds(i)  is  1  and the default value of valconds(i) 
%  is taken.
%  If  valconds  is not specified, then the default value for valconds(i) is
%       deriv. of cubic interpolant to nearest four points, if  conds(i)=1;
%       0                                                   if  conds(i)=2.

% C. de Boor / latest change: December 3, 1990
% C. de Boor / latest change: February 22, 1991 (added comments)
% C. de Boor / latest change: 12 April 1992 (added periodic case and comments)
% Copyright (c) 1990-92 by Carl de Boor and The MathWorks, Inc.


%     Generate the cubic spline interpolant in pp form.

if (nargin<3), conds=[1 1];end
if (nargin<4), valconds=[0 0]; end

n=length(x);[xi,ind]=sort(x);xi=xi(:);
if n<2,
   error('There should be at least two data points!')
elseif all(diff(xi))==0,
   error('The data abscissae should be distinct!')
elseif n~=length(y),
   error('Abscissa and ordinate vector should be of the same length!')
else   
   yi=y(ind);yi=yi(:);
      % set up the linear system for solving for the slopes at  xi .
   dx=diff(xi);divdif=diff(yi)./dx; 
   c = diag([dx(2:n-1);0],-1)+2*diag([0;dx(2:n-1)+dx(1:n-2);0]) ...
                 + diag([0;dx(1:n-2)],1);
   b=zeros(n,1);
   b(2:n-1)=3*(dx(2:n-1).*divdif(1:n-2)+dx(1:n-2).*divdif(2:n-1));
   if (~any(conds)),
      c(1,1)=1; c(1,n)=-1; b(1) = 0;
   elseif (conds(1)==2),
      c(1,1:2)=[2 1]; b(1)=3*divdif(1)-dx(1)*valconds(1)/2;
   else, 
      c(1,1:2) = [1 0]; b(1) = valconds(1);
      if (nargin<4|conds(1)~=1), % if endslope was not supplied, 
                                 % get it by local interpolation
         b(1)=divdif(1);
	 if (n>2), ddf=(divdif(2)-divdif(1))/(xi(3)-xi(1));
		   b(1) = b(1)-ddf*dx(1); end
	 if (n>3), ddf2=(divdif(3)-divdif(2))/(xi(4)-xi(2));
		   b(1)=b(1)+(ddf2-ddf)*(dx(1)*(xi(3)-xi(1)))/(xi(4)-xi(1));end
      end
   end
   if (~any(conds)),
      c(n,1:2)=dx(n-1)*[2 1]; c(n,n-1:n)= c(n,n-1:n)+dx(1)*[1 2];
      b(n) = 3*(dx(n-1)*divdif(1) + dx(1)*divdif(n-1));
   elseif (conds(2)==2),
      c(n,n-1:n)=[1 2]; b(n)=3*divdif(n-1)+dx(n-1)*valconds(2)/2;
   else, 
      c(n,n-1:n) = [0 1]; b(n) = valconds(2);
      if (nargin<4|conds(2)~=1), % if endslope was not supplied,
				 % get it by local interpolation
         b(n)=divdif(n-1);
	 if (n>2), ddf=(divdif(n-1)-divdif(n-2))/(xi(n)-xi(n-2));
	   b(n) = b(n)+ddf*dx(n-1); end
	 if (n>3), ddf2=(divdif(n-2)-divdif(n-3))/(xi(n-1)-xi(n-3));
	   b(n)=b(n)+(ddf-ddf2)*(dx(n-1)*(xi(n)-xi(n-2)))/(xi(n)-xi(n-3));end
      end
   end

     % solve for the slopes and convert to pp form
     % The final version should make use of MATLAB's eventual banded matrix
     % capability, by setting up a routine for complete cubic spline inter-
	 % polation to vector data, which can then be used, without pivoting,
	 % to solve for three splines, one matching function values and zero
	 % end slopes, while the other two match zero values and slopes except
	 % except for a slope of  1  at one endpoint or the other. Any side 
	 % condition can then be obtained as an appropriate linear combination
	 % of these three, with coefficients determined from a 2x2 linear 
	 % system.

   s=c\b;
   c4=(s(1:n-1)+s(2:n)-2*divdif(1:n-1))./dx;
   c3=(divdif(1:n-1)-s(1:n-1))./dx - c4;
   pp=ppmak_(xi',[c4./dx c3 s(1:n-1) yi(1:n-1)],1);
end
