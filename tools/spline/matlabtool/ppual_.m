function v=ppual_(pp,xx)
% PPUAL	Evaluate a pp.
%
%      v = ppual(pp,xx)
%
%  returns the value of the pp function  pp  at  xx.

% C. de Boor / latest change: June 19, 1989
% C. de Boor / latest change: December 1, 1990 (correct misprint in last line)
% Copyright (c) 1990-92 by Carl de Boor and The MathWorks, Inc.


%  if necessary, sort  xx 
xs=xx(:)';ix=length(xs);tosort=0;
if (length(find(diff(xs)<0))>0),
   tosort=1;[xs,indx]=sort(xs);
end

%  take apart  pp
[x,c,l,k,d]=ppbrk_(pp);

% for each data point, compute its breakpoint interval
index=max([sorted_(x(1:l),xs);ones(1,ix)]);

% now go to local coordinates
xs=xs-x(index);
% ... and apply nested multiplication:
if (d>1),
   ones(d,1)*xs; xs=ans(:)';
   1+d*ones(d,1)*index+[-d:-1]'*ones(1,ix); index=ans(:);
end
vv=c(index,1)';
for i=2:k,
   vv = xs.*vv + c(index,i)';
end

v=zeros(d,length(xx));v(:)=vv;
if tosort>0,[junk,indx]=sort(indx); v=v(:,indx);end

