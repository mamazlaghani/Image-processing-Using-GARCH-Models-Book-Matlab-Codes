function [LLF , G , H , e , h] = garchllfn_GG_noisy(Parameters , y , R , M , P , Q , van, X);
%GARCHLLFN Univariate GARCH process objective function (Gaussian innovations).
%   Compute the log-likelihood objective function value suitable for maximum 
%   likelihood estimation (MLE). This is the objective function optimized by 
%   the MATLAB Optimization Toolbox function FMINCON. The innovations are 
%   inferred from a conditional mean specification of ARMAX form, and the 
%   conditional variance of the innovations is then fit to a GARCH model. 
%   Gaussian innovations are assumed.
%
%   [LogLikelihood,Innovations,Sigma] = garchllfn(Parameters,Series,R,M,P,Q,X)
%
% Inputs:
%   Parameters - Column vector of process parameters associated with fitting
%     conditional mean and variance specifications to an observed return series 
%     Series. The conditional mean contributes the first (1 + R + M + nX) 
%     parameters (nX is the number of explanatory variables included in the 
%     regression component of the conditional mean); the conditional variance 
%     contributes the remaining (1 + P + Q) parameters. The length of Parameters 
%     is thus (2 + R + M + nX + P + Q).
%
%   Series - Matrix of observations of the underlying univariate return series 
%     of interest. Series is the response variable representing the time series 
%     fit to conditional mean and variance specifications. Each column of Series
%     is an independent realization (i.e., path). The last row of Series holds 
%     the most recent observation of each realization.
%
%   R - Non-negative, scalar integer representing the AR-process order.
%
%   M - Non-negative, scalar integer representing the MA-process order.
%
%   P - Non-negative, scalar integer representing the number of lags of the
%     conditional variance included in the GARCH process.
%
%   Q - Non-negative, scalar integer representing the number of lags of the 
%     squared innovations included in the GARCH process.
%
%   X - Time series regression matrix of explanatory variable(s). Typically, X 
%     is a regression matrix of asset returns (e.g., the return series of an 
%     equity index). Each column of X is an individual time series used as an 
%     explanatory variable in the regression component of the conditional mean. 
%     In each column of X, the first row contains the oldest observation and 
%     the last row the most recent. If empty or missing, the conditional mean 
%     will have no regression component
%
% Outputs:
%   LogLikelihood - Vector of log-likelihood objective function values evaluated
%     at the values in Parameters. The length of LogLikelihood is the same as 
%     the number of columns in Series. This function is meant to be optimized 
%     via the FMINCON function of the MATLAB Optimization Toolbox. FMINCON is 
%     a minimization routine, so LogLikelihood is actually the negative of what
%     is formally presented in most econometrics references.
%
%   G - Empty matrix (placeholder for forward compatibility).
%
%   H - Empty matrix (placeholder for forward compatibility).
%
%   Innovations - Innovations matrix inferred from the input Series matrix.
%
%   Sigma - Conditional standard deviation matrix corresponding to Innovations.
%
% Note: 
%   This is a performance-sensitive objective function called by an iterative
%   numerical optimizer. For this reason, no argument checking is performed.
%   Although this function may be called directly, it is better accessed by
%   calling GARCHINFER. Type "help garchinfer" for details.
%
% See also GARCHSIM, GARCHPRED, GARCHFIT, GARCHINFER.

% Copyright 1999-2002 The MathWorks, Inc.   
% $Revision: 1.9 $   $ Date: 1998/01/30 13:45:34 $

%
% References:
%   Bollerslev, T. (1986), "Generalized Autoregressive Conditional 
%     Heteroskedasticity", Journal of Econometrics, vol. 31, pp. 307-327.
%   Box, G.E.P., Jenkins, G.M., Reinsel, G.C., "Time Series Analysis: 
%     Forecasting and Control", 3rd edition, Prentice Hall, 1994.
%   Engle, Robert (1982), "Autoregressive Conditional Heteroskedasticity 
%     with Estimates of the Variance of United Kingdom Inflation", 
%     Econometrica, vol. 50, pp. 987-1007.
%   Hamilton, J.D., "Time Series Analysis", Princeton University Press, 1994.
%

%
% NOTES:
% Let y(t) = return series of interest (assumed stationary)
%     e(t) = innovations of the model noise process (assumed invertible)
%     h(t) = conditional variance of the innovations process e(t)
%
% The input coefficient vector 'Parameters' is formatted exactly as the 
% coefficients would be read from the recursive difference equations when 
% solving for the current values of the y(t) and h(t) time series. 
%
% Consider the following general ARMAX(R,M,nX)/GARCH(P,Q) form:
%
%   y(t) =  C + AR(1)y(t-1) + ... + AR(R)y(t-R) + e(t) 
%             + MA(1)e(t-1) + ... + MA(M)e(t-M) + B(1)X(t,1) + ... + B(nX)X(t,nX)
%
%   h(t) =  K + GARCH(1)h(t-1)   + ... + GARCH(P)h(t-P) 
%             +  ARCH(1)e(t-1)^2 + ... +  ARCH(Q)e(t-Q)^2
%
% For example, consider the following equations for the conditional mean 
% and variance of an ARMAX(R=2,M=2,nX=1)/GARCH(P=2,Q=2) composite model:
%
%   y(t) =  1.3 + 0.5y(t-1) - 0.8y(t-2) + e(t) - 0.6e(t-1) + 0.08e(t-2) + 1.2X(t)
%   h(t) =  0.5 + 0.2h(t-1) + 0.1h(t-2) + 0.3e(t-1)^2  +  0.2e(t-2)^2
%
% For the example listed above, the coefficient vector 'Parameters' would be
%
%   Parameters = [ C     AR(1:R)     MA(1:M)  B(1:nX)  K   GARCH(1:P)   ARCH(1:Q)]'
%              = [1.3   0.5 -0.8   -0.6 0.08    1.2   0.5    0.2 0.1     0.3 0.2 ]'
%
% Notice that the coefficient of e(t) in the conditional mean equation is
% defined to be 1, and is NOT included in 'Parameters' vector because it
% is not estimated.
%

%
% The mean coefficients are modified to accommodate the inference of the 
% innovations e(t) (i.e., the general conditional mean equation given 
% above for y(t) is now solved for e(t) as the dependent variable). The
% innovations are inferred from the following general conditional mean 
% equation:
%
%   e(t) = -C + y(t) - AR(1)y(t-1) - ... - AR(R)y(t-R)   
%                    - MA(1)e(t-1) - ... - MA(M)e(t-M) 
%                    -  B(1)X(t,1) - ... - B(nX)X(t,nX)
%
% In the code below, the ARMA(R,M) component of the conditional mean is
% de-coupled from the regression component for performance purposes. 
%
% Note that the input return series y(t) may have several columns in which
% each column represents a unique stochastic realization (i.e, each column
% is an independent path, or trial, of the process y(t)). In contrast, the 
% regression matrix X is NOT path-dependent in any way, and is simply a 
% known regression matrix applied to each realization of y(t).
%
% When extracting the ARMA component, coefficients are ordered as inferred 
% from the following ARMA component of the conditional mean:
%
%   e(t) = -C + y(t) - AR(1)y(t-1) - ... - AR(R)y(t-R)   
%                    - MA(1)e(t-1) - ... - MA(M)e(t-M) 
%

armaCoefficients  =  [-Parameters(1) ; 1 ; -Parameters(2:(1 + R + M))]';

%
% To infer the innovations e(t), we adopt the Box & Jenkins approach for
% conditioning. The iteration required to infer e(t) from the observed return
% series y(t) is started at date t = R + 1. The iteration is thus conditioned
% on the first R actual values of y(t), while the initial values of e(t) are 
% set to their expected value, E[e(t)] = 0. Note that, after the iteration 
% is complete, the first R elements of the innovations vector e(t) will be 0.
%

[T,nPaths]  =  size(y);          % T = # of samples, nPaths = # of realizations.

if isempty(X)          

%
%  General ARMA(R,M) form for conditional mean (no regression component).
%
   e  =  garchllfarmax(R , M , nPaths , armaCoefficients , y , T);

else                   

%
%  General ARMAX(R,M,nX) form for conditional mean.
%
%  Extract the coefficients of the regression component of the conditional 
%  mean. When extracting the regression component, coefficients are ordered
%  as inferred from the following regression component of the conditional mean:
%
%  e(t) = -B(1)X(t,1) - ... - B(nX)X(t,nX)
%

   regressCoefficients  =  -Parameters((2 + R + M):(1 + R + M + size(X,2)));

   e  =  garchllfarmax(R , M , nPaths , armaCoefficients , y , T , regressCoefficients , X);

end

%
% At this point, the innovations process e(t) has been inferred from the 
% observed return series y(t). By adopting the Box & Jenkins approach 
% outlined above, the first R values of e(t) are zero at this point. 
% Since GARCH modeling of the conditional variance requires the square
% of the innovations, e(t)^2, there is a subtle discrepancy between the
% conditional mean and variance processes.
%
% GARCH(P,Q) processing requires P pre-sample lags of h(t) and Q pre-sample
% lags of e(t)^2. Bollerslev suggests estimating the unconditional variance 
% of e(t) and assigning the estimate to the first max(P,Q) samples of h(t) 
% and e(t)^2. In the code below, a new vector, e2(t), is created whose
% initial (max(P,Q) + R) values are set to the unconditional variance of 
% e(t); the remaining values of e2 are just the squared innovations e(t)^2
% for t > R. Since the first R values of e(t) have been forced to zero, 
% the sample variance of e(t) is based on elements t > R (i.e., the 
% elements of e(t) that were actually computed, and NOT those that were 
% forced to zero as part of the conditioning) to help minimize transients.
% 

maxPQ    =  max(P,Q);
variance =  mean(e(R+1:end,:).^2);
e2       =  [variance(ones(maxPQ + R,1),:) ; e(R+1:end,:).^2];

%
% Update the number of samples to include the max(P,Q) pre-samples.
%

T  =  size(e2,1);

%
% Extract the conditional variance coefficients. 
%

varianceCoefficients  =  Parameters((2 + R + M + size(X,2)):end-2)';

%
% Over-write all conditional variance parameters less than or equal 
% to zero to prevent the LLF from becoming minus infinity or complex.
%

varianceCoefficients(varianceCoefficients <= 0) = realmin;

%
% Form the h(t) time series and evaluate the log-likelihood function.
% Note that LLF estimation ignores the pre-sample values of h(t) and
% e(t)^2 used for conditioning.
%

h    =  garchllfgarch(P , Q , nPaths , varianceCoefficients , e2 , T);
NNN=length(Parameters);
s=Parameters(NNN-1);
p=Parameters(NNN);

t    =  (maxPQ + 1):T;
e=e2.^.5;
%e=e((maxPQ + 1):T);
%h=h((maxPQ + 1):T);
LF = LF_fun_GG(e,van,h,s,p);
%LLF=0.5*sum(log(h(t,:))) + (T - maxPQ)*log([(s^2)*gamma(3/p)]/gamma(1/p)) +  sum(abs((e(t,:))./[(h(t,:).^.5)*s]).^p);
LLF=1;
%LLF=-(log(LF));
%LLF=1/LF;
%LLF  =  0.5*sum(log(h(t,:))) + (T - maxPQ)*log([(s^2)*gamma(3/p)]/gamma(1/p)) +  sum(abs((e(t,:))./[(h(t,:).^.5)*s]).^p);
%LLF  =  0.5 * (LLF  +  (T - maxPQ)*log(2*pi));

%
% Strip off the initial pre-sample values of h(t) and convert the 
% conditional variance to a standard deviation for output purposes.
%

h  =  sqrt(h((maxPQ + 1):T,:));

%
% Catch conditions that produce anomalous log-likelihood function values.
% Typically, what happens is that input parameter values will result in
% an unstable inverse filter, which produces LLF = inf. This, in turn, will
% not allow FMINCON to update the iteration. By setting the LLF to a large, 
% but finite, value, we can safeguard against anomalies and allows the 
% iteration to continue.
%

%if isnan(LLF) | isinf(LLF) | ~isreal(LLF)
 %  LLF  =  1e+20;
 %end

G = [];  % Placeholder for forward compatibility.
H = [];  % Placeholder for forward compatibility.
