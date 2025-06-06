function [LLF , G , H , e , h] = garchlfn0(Parameters , y , R , M , P , Q , X);

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
h    =  garchlfgarch(P , Q , varianceCoefficients , e2 , T);
t    =  (maxPQ + 1):T;
LLF  =  sum(log(h(t,:))) + sum((e2(t,:))./h(t,:));
LLF  =  0.5 * (LLF  +  (T - maxPQ)*log(2*pi));
h  =  sqrt(h((maxPQ + 1):T,:));
if isnan(LLF) | isinf(LLF) | ~isreal(LLF)
   LLF  =  1e+20;
end

G = [];  % Placeholder for forward compatibility.
H = [];  % Placeholder for forward compatibility.
