function y=diagonal(a);
N=length(a);
[b,d]=spdiags(a);
[n,m]=size(b);
y=zeros(N^2,1);
s=1;
for i=1:n*m
    if b(i)~=0
        y(s)=b(i);
        s=s+1;
    end
end
