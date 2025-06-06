function [y,yy]=weibull(a,b,n);
x=weibrnd(a,b,1,n);
[y,yy]=weibfit(x);
