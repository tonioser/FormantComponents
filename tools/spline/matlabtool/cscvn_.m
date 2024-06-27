function cs = cscvn_(points)
% CSCVN	Generate an interpolating parametric cubic spline curve.
%
%        cs = cscvn(points)
%
% returns a parametric `natural' cubic spline which interpolates to the given 
% points points(:,i)  at parameter values  t(i) , i=1,2,..., with  t(i)  chosen 
% by Eugene Lee's centripetal scheme, i.e., as accumulated squareroot of 
% chord-length.

% C. de Boor / latest change: Jan.28, 1990
%              latest change: May 12, 1991 change from csapn to csape
% Copyright (c) 1990-92 by Carl de Boor and The MathWorks, Inc.

[d,n] = size(points);

t=cumsum([0;((diff(points').^2)*ones(d,1)).^(1/4)])';

cs = csape_(t,points(1,:),[0 0]);

if (d>1),
   [breaks,coef,l,k] = ppbrk_(cs);
   coefs=zeros(d,l*k);
   coefs(1,:)=coef(:)';
   for dd=2:d
      [breaks,coef] = ppbrk_(csape_(t,points(dd,:),[0 0]));
      coefs(dd,:)=coef(:)';
   end
   temp = zeros(d*l,k); temp(:) = coefs;
   cs = ppmak_(breaks,temp,d);
end
