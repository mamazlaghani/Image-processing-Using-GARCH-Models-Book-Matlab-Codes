function [h,ee] = Egarchlfgarch2(p1 , p2 ,q1, q2 , nPaths, varianceCoefficients , e , variance, teta, gama,T1,T2,r1,r2,e2 );
maxpq1=max(p1,q1);
maxpq2=max(p2,q2);
g=zeros(T1,T2);

for i=1:maxpq2+r2
    g(:,i)=(variance);
end

for i=1:maxpq1+r1
    g(i,:)=variance;
end

h=zeros(T1,T2);
lnh=zeros(T1,T2);



for i=1:maxpq1
      lnh(i,:)=log(e2(i,:));   
end
for i=1:maxpq2
    lnh(:,i)=log(e2(:,i));
end
 %var=[0 varianceCoefficients(2:((p1+1)*(p2+1))) 0 varianceCoefficients(((p1+1)*(p2+1)+1):(p1+1)*(p2+1)+2) 1 varianceCoefficients(((p1+1)*(p2+1)+3):end)];
 var=[0 varianceCoefficients(2:((p1+1)*(p2+1))) 0 varianceCoefficients(((p1+1)*(p2+1)+1):end)];
 var1=var(1:((p1+1)*(p2+1)));
 var2=var(((p1+1)*(p2+1)+1):end);
 G=reshape(var1,p1+1,p2+1);
 A=reshape(var2,q1+1,q2+1);
 %%%%%
 %%%%%
 %
 %
 %
 %
 %
 for t=maxpq1+1:T1
     for k=maxpq2+1:T2
         %lnh(t,k)=varianceCoefficients(1)*[(1/(.01*exp(k)+.01*exp(t)))];
         %lnh(t,k)=varianceCoefficients(1)*[-(.1*t+.1*k)];
         %lnh(t,k)=varianceCoefficients(1)*[-(2.5*t+2.5*k)];
         lnh(t,k)=varianceCoefficients(1)*[(1/(2.5*exp(k)+2.5*exp(t)))];
         %lnh(t,k)=varianceCoefficients(1)+log(1+A1*t+A2*k);
         for i=1:p1+1
             for j=1:p2+1     
                 lnh(t,k)=lnh(t,k)+G(i,j)*lnh(t-i+1,k-j+1);
             end
         end
         for i=1:q1+1
             for j=1:q2+1                 
                 lnh(t,k)=lnh(t,k)+A(i,j)*g(t-i+1,k-j+1);
             end
         end
     end
    h(t,k)=exp(lnh(t,k));
    ee(t,k)=e(t-maxpq1,k-maxpq2)/(h(t,k)^.5);
    g(t,k)=teta*ee(t,k)+gama*(abs(ee(t,k))-0.7979);    
 end
 h=exp(lnh);
         