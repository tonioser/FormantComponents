function points=eval_sp(f,npoints)
% Returns npoints for the spline f
%     points=eval_sp(f,npoints);

if (f(1)==11), f=sp2pp_(f);end

[breaks,coefs,l,k,d]=ppbrk_(f);
x=breaks(1)+[0:npoints]*((breaks(l+1)-breaks(1))/npoints);
v=ppual_(f,x);

if (d==1),points=[x;v];
else,points=v([1,2],:);
end
points = points';