function [LLF , G , H , e , hfinal] = garchlfnGM1(Parameters , y , pp , R , M , P , Q , X, prob);

armaCoefficients  =  [-Parameters(1) ; 1 ; -Parameters(2:(1 + R + M))]';% parameters are arranged for innovation
[T,nPaths]  =  size(y); 
 e  =  garchlfarmax(R , M , nPaths , armaCoefficients , y , T);% e is innovation
 maxPQ    =  max(P,Q);
variance =  mean(e(R+1:end,:).^2);% initial value of conditional variance as described in main GARCH article

e2       =  [variance(ones(maxPQ + R,1),:) ; e(R+1:end,:).^2];
T  =  size(e2,1);
varianceCoefficients  =  Parameters((2 + R + M + size(X,2)):end)';
varianceCoefficients(varianceCoefficients <= 0) = realmin;
%h    =  garchllfgarch(P , Q , nPaths , varianceCoefficients , e2 , T);
for i=1:pp
h(:,i) = garchlfgarch0(P , Q ,nPaths, varianceCoefficients((P+Q+1)*(i-1)+1:(P+Q+1)*i) , e2 , T);
end
%prob=varianceCoefficients((P+Q+1)*pp+1:end);
%t    =  (maxPQ + 1):T;
%LLF  =  sum(log(h(t,:))) + sum((e2(t,:))./h(t,:));
%LLF  =  0.5 * (LLF  +  (T - maxPQ)*log(2*pi));
%LLF=0;
%for i=(maxPQ+1):T
%    LLF=LLF+.5*(log(2*pi)+log(h(i))+e2(i)/h(i));
%end 
e1=sqrt(e2);
h;
LLF=0;
for i=(maxPQ+1):T
   LLF=LLF-log(gaussianmixture(e1(i),prob,zeros(pp,1),h(i,:)));
end 
[aaa,bbb]=size(h);
sx=zeros(aaa,1);
for i=1:pp
    sx=sx+(prob(i)*h(:,i));
end
hhfinal=sx;
hfinal=sqrt(hhfinal((maxPQ + 1):T,:));
%h  =  sqrt(h((maxPQ + 1):T,:));
%if isnan(LLF) | isinf(LLF) | ~isreal(LLF)
 
%LLF  =  1e+20;
%end

G = [];  % Placeholder for forward compatibility.
H = [];  % Placeholder for forward compatibility.
