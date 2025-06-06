function xx=generalized_laplacian_rnd(s,p,N1,N2,lowerlimit,upperlimit,M);
%lowerlimit and upperlimit are related to inverse cdf of GG and are the
%minx and maxx 
%N is the number of RVs
%M is the number uses in cdf_inverse (100) for better precision this number
%must be large
n=rand(N1,N2);
xx=zeros(N1,N2);
for i=1:N1
    for j=1:N2
        xx(i,j)=generalized_laplacian_cdf_inverse(n(i,j),s,p,M,lowerlimit,upperlimit);    
    end
end
