function h=garchlfgarch(P , Q ,npaths, varianceCoefficients , e2 , T);
for i=1:max(P,Q)
    h(i)=e2(i);
end
ee2=e2(max(P,Q)+1:end);
%size(ee2)
%T
for i=max(P,Q)+1:T-P+1
                  h(i)= varianceCoefficients(1)+varianceCoefficients(2:P+1)*h(i-P:i-1)'+varianceCoefficients(P+2:end)*ee2(i-Q:i-1);
end
for i=T-P+1:T
    h(i)=h(T-P+1);
end
   h=h'; 