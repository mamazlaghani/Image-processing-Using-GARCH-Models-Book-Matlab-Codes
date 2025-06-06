function x=lognormal(m,s,n);
data=lognrnd(m,s,1,n);
x=fminsearch('lognormalmin',[0 1],'',data);
