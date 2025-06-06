function x=normal(m,n);
data=normrnd(m,1,1,n);
x=fminsearch('normalmin',0,'',data);
