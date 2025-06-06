function e=garchlfarmax28(r1, r2 , m1 , m2 , nPaths , armaCoefficients , y , T1, T2);
e=zeros(T1+m1,T2+m2);
%
%
%
arma=[0 armaCoefficients(3:([(r1+1)*(r2+1)]+1)) 0 armaCoefficients(((r1+1)*(r2+1)+2):end)];
size(arma);
arma1=arma(1:((r1+1)*(r2+1)));
 arma2=arma(((r1+1)*(r2+1)+1):end);
 size(arma2);
 m1+1;
 m2+1;
 AR=reshape(arma1,r1+1,r2+1);
 MA=reshape(arma2,m1+1,m2+1);
 %
 %
 %
 for i=1+m1:m1+r1
     for j=1+m2:m2+r2
         e(i,j)=0;
     end
 end
 for t=m1+r1+1:m1+T1
     for k=m2+r2+1:m2+T2
         e(t,k)=armaCoefficients(1)+y(t-m1,k-m2);
         for i=1:r1+1
             for j=1:r2+1
                 e(t,k)=e(t,k)+AR(i,j)*y(t-i+1-m1,k-j-m2+1);
             end
         end
         for i=1:m1+1
             for j=1:m2+1
                 e(t,k)=e(t,k)+MA(i,j)*e(t-i+1,k-j+1);
             end
         end
     end
 end
 e=e(m1+1:end,m2+1:end);
 