function [LLF , G , H , e , hh] = garchlfn22_GM(Parameters , y , r1 , r2 , m1 , m2 , p1 ,p2 ,q1 , q2 , X,pnum);

%armaCoefficients  =  [-Parameters(1) ; 1 ; -Parameters(2:(1 + r1 + m1 + r2 + m2))]';% parameters are arranged for innovation
armaCoefficients  =  [-Parameters(1) ; 1 ; -Parameters(2:(1 + r1)*(r2+1)-1+(m1 + 1)*(m2+1)-1+1)]';% parameters are arranged for innovation
[T1,T2]=size(y);
%yy=reshape(y,a*b,1);
%[T,nPaths]  =  size(yy); 
nPaths=1;
 e  =  garchlfarmax2(r1 , r2 , m1 , m2 , nPaths , armaCoefficients , y , T1 , T2 );
 %e=reshape(ee,a,b);
 maxpq1=max(p1,q1);
 maxpq2=max(p2,q2);
 %variance = mean( mean(e(r1+1:end,:).^2));% initial value of conditional variance as described in main GARCH article
  variance = mean( mean(e(r1+1:end,r2+1:end).^2));% initial value of conditional variance as described in main GARCH article
   
 for i=1:maxpq2+r2
      e2(:,i)=variance;
 end
 for i=1:maxpq1+r1
     e2(i,:)=variance;
 end
 for i=maxpq1+1+r1:T1+maxpq1
     for j=maxpq2+1+r2:T2+maxpq2
         e2(i,j)=e(i-maxpq1,j-maxpq2).^2;
     end
 end
[T1,T2]  =  size(e2);
varianceCoefficients  =  Parameters((2 + r1 + r2 + m1 + m2 + size(X,2)):end-(pnum*2))';
varianceCoefficients(varianceCoefficients <= 0) = realmin;
h    =  garchlfgarch2(p1 , p2 ,q1, q2 , nPaths, varianceCoefficients , e2 , T1 , T2);
prob=Parameters(end-2*pnum+1:end-pnum);
mean1=Parameters(end-pnum+1:end);
e1=sqrt(e2);
LLF=0;
for i=maxpq1+1:T1
    for j=maxpq2+1:T2
         LLF=LLF-log(gaussianmixture(e1(i,j),prob,mean1,h(i,j)*ones(1,pnum)));
     end
 end
 
hh=sqrt(h(maxpq1+1:T1,maxpq2+1:T2)); 
 if isnan(LLF) | isinf(LLF) | ~isreal(LLF)
   LLF  =  1e+20;
end

G = [];  % Placeholder for forward compatibility.
H = [];  % Placeholder for forward compatibility.
