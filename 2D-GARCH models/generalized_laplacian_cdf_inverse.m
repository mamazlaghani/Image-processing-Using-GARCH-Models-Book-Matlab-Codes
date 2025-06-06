function xx=generalized_laplacian_cdf_inverse(f,s,p,N,lowerlimit,upperlimit);
x=linspace(lowerlimit,upperlimit,N);
for i=1:N
    y(i)=generalized_laplacian_cdf(x(i),s,p);
end
for i=1:N
    if (y(i)>f) | (y(i)==f)
        xx=(x(i)+x(i-1))/2;
        break;
    end
end