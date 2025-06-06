function [a,b]=beta(x,y,n);
r=betarnd(x,y,1,n);
[a,b]=betafit(r);