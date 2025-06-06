function f = dft2(x,w1,w2);
s1=0;
[n,m]=size(x);
for t=1:n
    for k=1:m
        s1 = s1+x(t,k)*(exp(-i*w1*t))*(exp(-i*w2*k));
    end
end
f=s1;