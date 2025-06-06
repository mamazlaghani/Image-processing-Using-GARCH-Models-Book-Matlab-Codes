function x=fminn(b);
c=b;
f=inline('(x-c)^2');
x=fminsearch(f,[1,1]);