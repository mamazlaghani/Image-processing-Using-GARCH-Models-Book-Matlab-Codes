function [xx,h]=wienerin(s,x);
%s=normrnd(0,1,1,n);
cx=xcorr(x);
cxs=xcorr(x,s);
a=toeplitz(cx,cx);
b=cxs';
h=lsqr(a,b);
xx=conv(h,x);
