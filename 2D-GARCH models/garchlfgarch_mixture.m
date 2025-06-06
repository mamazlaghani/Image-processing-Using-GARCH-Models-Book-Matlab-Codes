function [h,h1,h2]=garchlfgarch_mixture(P,Q, nPaths, varianceCoefficients1 , varianceCoefficients2, e2 , T, prob, e, m);
prob;
maxpq=max(P,Q);
h=zeros(T,1);
for i=1:maxpq
      h(i)=e2(i);   
end


%first element of mixture
h1=zeros(T,1);
for i=1:maxpq
      h1(i)=e2(i);   
end
%second element of mixture
h2=zeros(T,1);
for i=1:maxpq
      h2(i)=e2(i);   
end


%first element of mixture coefficients
 G1=varianceCoefficients1(2);
 A1=varianceCoefficients1(3);
 
 %second element of mixture coefficients
 G2=varianceCoefficients2(2);
 A2=varianceCoefficients2(2);
 
 %%%%%
 %%%%%
 %
 %
 %
 %
 %
 b=zeros(T,1);
 for t=maxpq+1:T
     h1(t)=varianceCoefficients1(1);
     h2(t)=varianceCoefficients2(1);
     for i=1:P
         h1(t)=h1(t)+G1(i)*h(t-i);
         h2(t)=h2(t)+G2(i)*h(t-i);
     end
     for i=1:Q
         h1(t)=h1(t)+A1(i)*e2(t-i);
         h2(t)=h2(t)+A2(i)*e2(t-i);
     end
        b(t)=prob*normpdf(e(t),m,h1(t)^.5)+(1-prob)*normpdf(e(t),m,h2(t)^.5);
     if b(t) == 0
         h(t)=prob*h1(t)+(1-prob)*h2(t);
     else
         h(t)=[(prob*normpdf(e(t),m,h1(t)^.5))/(prob*normpdf(e(t),m,h1(t)^.5)+(1-prob)*normpdf(e(t),m,h2(t)^.5))]*h1(t)+[((1-prob)*normpdf(e(t),m,h2(t)^.5))/(prob*normpdf(e(t),m,h1(t)^.5)+(1-prob)*normpdf(e(t),m,h2(t)^.5))]*h2(t);
     end
 end
         
 
         