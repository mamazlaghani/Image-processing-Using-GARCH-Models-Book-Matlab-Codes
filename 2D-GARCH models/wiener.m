function [xx,h]=wiener(x,s);
cx=xcorr(x,'coeff');
cxs=xcorr(x,s,'coeff');
a=toeplitz(cx);
b=cxs';
h=cgs(a,b);
xx=filter(h,1,x)*2;
%xx=conv(h,x);
