function [s,p]=generalized_laplacian(c);
r=c(:);
va=var(r);
k=[mean(r.^4)/(va^2)];
p=fminsearch('differnce',.6,[],k);
s2=(va*gamma(1/p))/(gamma(3/p));
s=s2^.5;