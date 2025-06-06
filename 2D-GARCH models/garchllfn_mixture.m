function [LLF , G , H , e , hh, hh1 , hh2]  = garchllfn_mixture(Parameters , y , R , M , P , Q , X)
armaCoefficients  =  [-Parameters(1) ; 1 ; -Parameters(2:(1 + R + M))]';
varianceCoefficients1  =  Parameters((2 + R + M + size(X,2)):(2 + R + M + size(X,2))+P+Q)';
varianceCoefficients2  =  Parameters((2 + R + M + size(X,2))+P+Q+1:(2 + R + M + size(X,2))+2*P+2*Q+1)';
prob=Parameters((2 + R + M + size(X,2))+2*P+2*Q+2);
[T,nPaths]  =  size(y);          % T = # of samples, nPaths = # of realizations.

e  =  garchllfarmax(R , M , nPaths , armaCoefficients , y , T);

maxPQ    =  max(P,Q);
variance =  mean(e(R+1:end,:).^2);
e2       =  [variance(ones(maxPQ + R,1),:) ; e(R+1:end,:).^2];
yy       =  [variance(ones(maxPQ + R,1),:).^.5 ; y(R+1:end,:).^2];
T  =  size(e2,1);
varianceCoefficients1(varianceCoefficients1 <= 0) = realmin;
varianceCoefficients2(varianceCoefficients2 <= 0) = realmin;
%
% Form the h(t) time series and evaluate the log-likelihood function.
% Note that LLF estimation ignores the pre-sample values of h(t) and
% e(t)^2 used for conditioning.
%
mmm=Parameters(1);
prob;
[h,h1,h2]    =  garchlfgarch_mixture(P,Q , nPaths, varianceCoefficients1 , varianceCoefficients2, e2 , T, prob, yy, mmm);
prob;
LLF=0;
for i=maxPQ+1:T
    LLF=LLF-log([prob*[(1/((2*pi*h1(i))^.5)*exp(-e2(i)/(2*h1(i))))]+(1-prob)*[(1/((2*pi*h2(i))^.5)*exp(-e2(i)/(2*h2(i))))]]+.00001);
end

 
%
% Strip off the initial pre-sample values of h(t) and convert the 
% conditional variance to a standard deviation for output purposes.
%

hh  =  sqrt(h((maxPQ + 1):T,:));
hh1  =  sqrt(h1((maxPQ + 1):T,:));
hh2  =  sqrt(h2((maxPQ + 1):T,:));

%
% Catch conditions that produce anomalous log-likelihood function values.
% Typically, what happens is that input parameter values will result in
% an unstable inverse filter, which produces LLF = inf. This, in turn, will
% not allow FMINCON to update the iteration. By setting the LLF to a large, 
% but finite, value, we can safeguard against anomalies and allows the 
% iteration to continue.
%

if isnan(LLF) | isinf(LLF) | ~isreal(LLF)
   LLF  =  1e+20;
end

G = [];  % Placeholder for forward compatibility.
H = [];  % Placeholder for forward compatibility.
