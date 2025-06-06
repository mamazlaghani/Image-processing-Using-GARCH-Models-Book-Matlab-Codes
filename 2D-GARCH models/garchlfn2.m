function [LLF , G , H , e , hh] = garchlfn2(Parameters , y , r1 , r2 , m1 , m2 , p1 ,p2 ,q1 , q2 , X);

armaCoefficients  =  [-Parameters(1) ; 1 ; -Parameters(2:(1 + r1 + m1 + r2 + m2))]';% parameters are arranged for innovation
[a,b]=size(y);
yy=reshape(y,a*b,1);
[T,nPaths]  =  size(yy); 
nPaths=1;
 ee  =  garchllfarmax(r1 , m1 , nPaths , armaCoefficients , y , T );
 e=reshape(ee,a,b);
 maxpq1=max(p1,q1);
 maxpq2=max(p2,q2);
 variance = mean( mean(e(r1+1:end,:).^2));% initial value of conditional variance as described in main GARCH article
   for i=1:maxpq2
      e2(:,i)=variance;
 end
 for i=1:maxpq1
     e2(i,:)=variance;
 end
 for i=maxpq1+1:a+maxpq1
     for j=maxpq2+1:b+maxpq2
         e2(i,j)=e(i-maxpq1,j-maxpq2).^2;
     end
 end

[T1,T2]  =  size(e2);
varianceCoefficients  =  Parameters((2 + r1 + r2 + m1 + m2 + size(X,2)):end)';
varianceCoefficients(varianceCoefficients <= 0) = realmin;
h    =  garchlfgarch2(p1 , p2 ,q1, q2 , nPaths, varianceCoefficients , e2 , T1 , T2);
LLF=0;
for i=maxpq1+1:T1
    for j=maxpq2+1:T2
         LLF=LLF+.5*(log(2*pi)+log(h(i,j))+e2(i,j)/h(i,j));
     end
 end
 
hh=sqrt(h(maxpq1+1:T1,maxpq2+1:T2)); 
 if isnan(LLF) | isinf(LLF) | ~isreal(LLF)
   LLF  =  1e+20;
end

G = [];  % Placeholder for forward compatibility.
H = [];  % Placeholder for forward compatibility.
