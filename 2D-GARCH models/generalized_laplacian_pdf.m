function f=generalized_laplacian_pdf(x,s,p);
z=2*(s/p)*gamma(1/p);
N=length(x);
f=zeros(1,N);
for i=1:N
    f(i)=[exp(-1*abs((x(i)/s)^p))]/z;
end
