function y=dft(x,n,f,k)
s=0;
for j=1:n
    s=x(j)*exp(-1*(j-1-k)*f*i)+s;
end
y=s;

    