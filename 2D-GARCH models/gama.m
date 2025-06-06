function [y,yy]=gama(a,b,n);
x=gamrnd(a,b,1,n);
[y,yy]=gamfit(x);
