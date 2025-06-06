function [K , GARCH , ARCH] = garch0(P , Q , unconditionalVariance)
%GARCH0 Initial GARCH process parameter estimates.
%   Given the orders of a GARCH(P,Q) model and an estimate of the unconditional 
%   variance of the innovations process, compute initial estimates for the 
%   (1 + P + Q) parameters of a GARCH(P,Q) conditional variance model. These
%   estimates serve as initial guesses for further refinement via maximum 
%   likelihood.
%
%   [K , GARCH , ARCH] = garch0(P , Q , Variance)
%
% Inputs:
%   P - Non-negative, scalar integer representing the number of lags of the
%     conditional variance included in the GARCH process.
%
%   Q - Non-negative, scalar integer representing the number of lags of the 
%     squared innovations included in the GARCH process.
%
%   Variance - Estimate of the unconditional variance of the innovations noise
%     process. Equivalently, it may also be viewed as the variance estimate of
%     the white noise innovations under the assumption of homoskedasticity.
%
% Outputs:
%   K - Conditional variance constant (scalar).
%
%   GARCH - P-element column vector of coefficients of lagged conditional 
%     variances.
%
%   ARCH - Q-element column vector of coefficients of lagged squared 
%     innovations.
%
% Note: 
%   This is an internal helper function for GARCHFIT. No error checking is 
%   performed.
%

%
% The following initial guesses are based on empirical observation of GARCH
% model parameters. This approach is rather ad hoc, but is very typical of
% GARCH models in financial time series. The most common (and very useful!) 
% GARCH model is the simple GARCH(1,1) model in which the coefficient of the
% lagged conditional variance (i.e., the 'GARCH' coefficient) is about 0.8 
% to 0.9, and the coefficient of lagged squared innovation (i.e., the 'ARCH' 
% coefficient) is about 0.05. Thus, a reasonable GARCH(1,1) assumption is:
%
%            h(t) = K + 0.85h(t-1) + 0.05e^2(t-1)
%
% In a GARCH(1,1) model, the unconditional variance of the innovations 
% process, V, is
%
%            V = K / (1 - (0.85 + 0.05))
% or,
%
%            K = V*(1 - (0.85 + 0.05))
%
% For higher-order GARCH(P,Q) models, this approach assumes the sum of the
% lagged conditional variance coefficient = 0.85, and the sum of the coefficients
% of lagged squared innovations = 0.05.
%

GARCH =  0.85;
GARCH =  GARCH(ones(P,1)) / max(P,1);
GARCH =  GARCH(:);

ARCH  =  0.05;
ARCH  =  ARCH(ones(Q,1)) / max(Q,1);
ARCH  =  ARCH(:);


if isempty(unconditionalVariance) |  (unconditionalVariance <= 0)
   K  =  1e-5;   % A decent assumption for daily returns.
else
   K  =  unconditionalVariance * (1 - sum([GARCH ; ARCH]));
end