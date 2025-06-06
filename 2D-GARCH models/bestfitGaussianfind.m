function p=bestfitGaussianfind(r);
%N=length(r);
[h,x]=hist(r,200);
hh=h./(sum(h*(x(2)-x(1))));
m=mean(r);
st=std2(r);
p0=[m st];
p=fminsearch('bestfitGaussian'  , p0 ,[],x,hh);