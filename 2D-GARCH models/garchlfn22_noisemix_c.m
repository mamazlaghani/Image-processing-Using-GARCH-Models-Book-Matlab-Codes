function [LLF , G , H , e , hh, hh1 , hh2] = garchlfn22_noisemix_c(Parameters , y , r1 , r2 , m1 , m2 , p1 ,p2 ,q1 , q2 ,me1,me2,va1,va2,prob, X);

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
      yy(:,i)=variance^.5;
 end
 for i=1:maxpq1+r1
     e2(i,:)=variance;
     yy(i,:)=variance^.5;
 end
 for i=maxpq1+1+r1:T1+maxpq1
     for j=maxpq2+1+r2:T2+maxpq2
         e2(i,j)=e(i-maxpq1,j-maxpq2).^2;
         yy(i,j)=y(i-maxpq1,j-maxpq2);
     end
 end
 size(e);
 size(yy);
[T1,T2]  =  size(e2);
%T1%
%T2%
varianceCoefficients  =  Parameters((2 + r1 + r2 + m1 + m2 + size(X,2)):end)';
varianceCoefficients(varianceCoefficients <= 0) = realmin;
% varianceCoefficients  =  Parameters((2 + r1 + r2 + m1 + m2 + size(X,2)):(2 + r1 + r2 + m1 + m2 + size(X,2))+(p1+1)*(p2+1)+(q1+1)*(q2+1)-2)';
%varianceCoefficients(varianceCoefficients <= 0) = realmin;
mmm=Parameters(1);
[h,h1,h2]    =  garchlfgarch2_mixture_c(p1 , p2 ,q1, q2 , nPaths, varianceCoefficients , e2 , T1 , T2, prob, yy, mmm, va1,va2,me1,me2);

LLF=0;
e21=(yy-me1).^2;
e22=(yy-me2).^2;
for i=maxpq1+1:T1
    for j=maxpq2+1:T2
        %         LF=LF*([(1/((h1(i,j))^.5)*exp(-e2(i,j)/(2*h1(i,j))))]);
         %LLF=LLF+.5*(log(2*pi)+log(h(i,j))+e2(i,j)/h(i,j));
         LLF=LLF-log([prob*[(1/((2*pi*h1(i,j))^.5)*exp(-e21(i,j)/(2*h1(i,j))))]+(1-prob)*[(1/((2*pi*h2(i,j))^.5)*exp(-e22(i,j)/(2*h2(i,j))))]]+.00001);

     end
 end
% LLF=-log(LF);
%T1
%maxpq1
hh=sqrt(h(maxpq1+1:T1,maxpq2+1:T2)); 
hh1=sqrt(h1(maxpq1+1:T1,maxpq2+1:T2)); 
hh2=sqrt(h2(maxpq1+1:T1,maxpq2+1:T2)); 
hh(hh<= 0) = realmin;
hh1(hh1<= 0) = realmin;
hh2(hh2<= 0) = realmin;
 if isnan(LLF) | isinf(LLF) | ~isreal(LLF)
   LLF  =  1e+20;
end

G = [];  % Placeholder for forward compatibility.
H = [];  % Placeholder for forward compatibility.
