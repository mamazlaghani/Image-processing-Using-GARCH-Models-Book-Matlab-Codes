function f=generalized_laplacian_cdf(x,s,p);
N=length(x);
f=zeros(1,N);
T=zeros(1,N);
for i=1:N
    f(i)=quad('generalized_laplacian_pdf',-1000,x(i),1.e-6,[],s,p);
end
