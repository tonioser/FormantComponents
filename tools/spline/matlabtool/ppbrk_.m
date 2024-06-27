function [breaks,coefs,l,k,d]=ppbrk_(pp,print)
% PPBRK	Break a pp into pieces.
% 
%        [breaks,coefs,l,k,d] = ppbrk(pp[,print])
%
%  breaks apart the  pp  function into its pieces, and, optionally, prints them.

% C. de Boor / latest change: Feb.25, 1989
% Copyright (c) 1990-92 by Carl de Boor and The MathWorks, Inc.


if (pp(1)==10),
   d=pp(2);
   l=pp(3);
   breaks=pp(3+[1:l+1]);
   k=pp(5+l);
   coefs=zeros(d*l,k);coefs(:)=pp(5+l+[1:d*l*k]);
   if (nargin>1),breaks,coefs,l,k,d,end
else,
   error('the input array does not seem to describe a pp function')
end
