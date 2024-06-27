function pp=ppmak_(breaks,coefs,d)
% PPMAK	Make a pp.
%
%     pp = ppmak(breaks,coefs[,d])
%
%  puts together a pp function from the breaks and coefficients input or
%  requested. 
% 
%    If  d  is not explicitly given, program expects  coefs(d,(k,l)) ,
%  and infers  l := length(breaks)-1  and then  k := (# cols(coefs))/l .
% 
%    If  d  is explicitly given, program expects  coefs((d,l),k) , and
%  infers  k := # cols(coefs), and  l := (# rows(coefs))/d .

% C. de Boor / latest change: June 1, 1989
% Copyright (c) 1990-92 by Carl de Boor and The MathWorks, Inc.


if (nargin==0);
   breaks=input('Give the (l+1)-vector of breaks  >');  
   coefs=input('Give the (d by (k*l)) matrix of local pol. coefficients  >');
end
if (nargin<3),
   [d,kl]=size(coefs);
   if (d==0),
      error('The coefficient sequence is empty!')
   end
   l=length(breaks)-1;k=fix(kl/l);
   if (k<=0)|(k*l~=kl);
      fprintf('The given number %.0f of polynomial pieces is incompatible',l)
      fprintf(' with the total number %.0f of coefficients supplied!\n',kl)
      error('')
   elseif (~isempty(find(diff(breaks)<0))),
      error('The breakpoint sequence should be nondecreasing!')
   else
      % the pp-format expects coefs in array  (d*l) by k, while the standard
      % input supplies them in an array d by (k*l) . This requires the
      % following shuffling, from  D+d(-1+K + k(-1+L))=D-d +(K-k)d + dkL
      % to  D+d(-1+L + l(-1+K)=D-d +(L-l)d + dlK .
      c=coefs(:);([1-k:0]'*ones(1,l)+k*ones(k,1)*[1:l])';
      coefs=[1-d:0]'*ones(1,kl)+d*ones(d,1)*(ans(:)');
      coefs(:)=c(coefs);
   end
else
   [dl,k]=size(coefs); l = dl/d;
   if (l+1~=length(breaks)|dl~=l*d), 
      error('input arguments have incompatible sizes.'),return,end
end
pp=[10 d l breaks(:)' k coefs(:)'];
