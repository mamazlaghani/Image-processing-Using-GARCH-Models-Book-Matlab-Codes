function h=garchlfgarch0(P , Q ,npaths, varianceCoefficients , e2 , T);
for i=1:max(P,Q)
    h(i)=e2(i);
end
for i=max(P,Q)+1:T
    h(i)=varianceCoefficients(1);
    for j=1:P
        h(i)=h(i)+varianceCoefficients(1+j)*h(i-j);
    end
    for j=1:Q
        h(i)=h(i)+varianceCoefficients(1+P+j)*e2(i-j);
    end
end
h=h';
