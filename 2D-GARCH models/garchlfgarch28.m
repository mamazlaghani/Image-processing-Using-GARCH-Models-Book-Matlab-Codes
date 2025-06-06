function h=garchlfgarch28(varianceCoefficients, e2 , T1 , T2 );
maxpq1=max(p1,q1);
maxpq2=max(p2,q2);
h=zeros(T1,T2);
for i=1:maxpq1
      h(i,:)=e2(i,:);   
end
for i=1:maxpq2
    h(:,i)=e2(:,i);
end
 var=[0 varianceCoefficients(2:9) 0 varianceCoefficients(5:7)];
 a0=varianceCoefficients(1);
 a11=varianceCoefficients(2);
 a10=varianceCoefficients(3);
 a01=varianceCoefficients(4);
 am1m1=varianceCoefficients(5);
 am10=varianceCoefficients(6);
 a0m1=varianceCoefficients(7);
 am11=varianceCoefficients(8);
 a1m1=varianceCoefficients(9);
 
 
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
         h(t,k)=varianceCoefficients(1);
         for i=1:p1+1
             for j=1:p2+1     
                 h(t,k)=h(t,k)+G(i,j)*h(t-i+1,k-j+1);
             end
         end
         for i=1:q1+1
             for j=1:q2+1                 
                 h(t,k)=h(t,k)+A(i,j)*e2(t-i+1,k-j+1);
             end
         end
     end
 end
 
         