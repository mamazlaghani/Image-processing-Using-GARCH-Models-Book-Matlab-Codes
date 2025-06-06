function h=linear(y,a);
varianceCoefficients=a(2:end);
[T1,T2]=size(y);
e=y-a(1);
variance = mean( mean(e(1:end,1:end).^2));
for i=1:1
      e2(:,i)=variance;
 end
 for i=1:1
     e2(i,:)=variance;
 end
 for i=1+1:T1+1
     for j=1+1:T2+1
         e2(i,j)=e(i-1,j-1).^2;
     end
 end
[T1,T2]  =  size(e2);
varianceCoefficients(varianceCoefficients <= 0) = realmin;
nPaths=1;
h    =  garchlfgarch2(1 , 1 ,1, 1 , nPaths, varianceCoefficients , e2 , T1 , T2);
h=h(1+1:T1,1+1:T2)