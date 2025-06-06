function [h,h1,h2]=garchlfgarch2_mixture(p1 , p2 ,q1, q2 , nPaths, varianceCoefficients1 , varianceCoefficients2, e2 , T1 , T2, prob, e, m);
p=prob;
maxpq1=max(p1,q1);
maxpq2=max(p2,q2);
h=zeros(T1,T2);
for i=1:maxpq1
      h(i,:)=e2(i,:);   
end
for i=1:maxpq2
    h(:,i)=e2(:,i);
end

%first element of mixture
h1=zeros(T1,T2);
for i=1:maxpq1
      h1(i,:)=e2(i,:);   
end
for i=1:maxpq2
    h1(:,i)=e2(:,i);
end

%second element of mixture
h2=zeros(T1,T2);
for i=1:maxpq1
    h2(i,:)=e2(i,:);   
end
for i=1:maxpq2
    h2(:,i)=e2(:,i);
end


%first element of mixture coefficients
 varr=[0 varianceCoefficients1(2:((p1+1)*(p2+1))) 0 varianceCoefficients1(((p1+1)*(p2+1)+1):end)];
 var1=varr(1:((p1+1)*(p2+1)));
 var2=varr(((p1+1)*(p2+1)+1):end);
 G1=reshape(var1,p1+1,p2+1);
 A1=reshape(var2,q1+1,q2+1);
 
 %second element of mixture coefficients
 varr=[0 varianceCoefficients2(2:((p1+1)*(p2+1))) 0 varianceCoefficients2(((p1+1)*(p2+1)+1):end)];
 var1=varr(1:((p1+1)*(p2+1)));
 var2=varr(((p1+1)*(p2+1)+1):end);
 G2=reshape(var1,p1+1,p2+1);
 A2=reshape(var2,q1+1,q2+1);
 
 %%%%%
 %%%%%
 %
 %
 %
 %
 %
 b=zeros(T1,T2);
 for t=maxpq1+1:T1
     for k=maxpq2+1:T2
         h1(t,k)=varianceCoefficients1(1);
         h2(t,k)=varianceCoefficients2(1);
         for i=1:p1+1
             for j=1:p2+1     
                 h1(t,k)=h1(t,k)+G1(i,j)*h(t-i+1,k-j+1);
                 h2(t,k)=h2(t,k)+G2(i,j)*h(t-i+1,k-j+1);
             end
         end
         for i=1:q1+1
             for j=1:q2+1                 
                 h1(t,k)=h1(t,k)+A1(i,j)*e2(t-i+1,k-j+1);
                 h2(t,k)=h2(t,k)+A2(i,j)*e2(t-i+1,k-j+1);
             end
         end
         b(t,k)=prob*normpdf(e(t,k),m,h1(t,k)^.5)+(1-prob)*normpdf(e(t,k),m,h2(t,k)^.5);
         if b(t,k) == 0
            h(t,k)=prob*h1(t,k)+(1-prob)*h2(t,k);
         else
            h(t,k)=[(prob*normpdf(e(t,k),m,h1(t,k)^.5))/(prob*normpdf(e(t,k),m,h1(t,k)^.5)+(1-prob)*normpdf(e(t,k),m,h2(t,k)^.5))]*h1(t,k)+[((1-prob)*normpdf(e(t,k),m,h2(t,k)^.5))/(prob*normpdf(e(t,k),m,h1(t,k)^.5)+(1-prob)*normpdf(e(t,k),m,h2(t,k)^.5))]*h2(t,k);
         end

     end
 end
 
         