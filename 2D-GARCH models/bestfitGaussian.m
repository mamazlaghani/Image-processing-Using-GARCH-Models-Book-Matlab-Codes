function r=bestfitGaussian(p,x,h);
%N=length(h);
m=p(1);
st=p(2);
f=normpdf(x,m,st);
r=sum((f-h).^2);


