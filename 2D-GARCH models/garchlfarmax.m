 function e=garchlfarmax(R , M , nPaths , armaCoefficients , y , T);
 e=zeros(T+M,1);%set all values to zero
 %set M initial value of e to zero
 for i=1+M:M+R
     e(i)=0;
 end
 %calculating innovation through the ARMA formulation
 for i=M+R+1:M+T
     e(i)=armaCoefficients(1)+y(i-M);
     for j=1:R
         e(i)=e(i)+armaCoefficients(j+2)*y(i-j-M);
     end
     for j=1:M
         e(i)=e(i)+armaCoefficients(j+2+R)*e(i-j);
     end
 end
 e=e(M+1:end);
     
     