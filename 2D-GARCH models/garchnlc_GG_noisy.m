function [c, ceq, gc, gceq] = garchnlc(Parameters , y , R , M , P , Q ,van, X)
%GARCHNLC GARCH Toolbox non-linear constraints function.
%   Enforce any non-linear constraints needed for GARCH maximum likelihood
%   parameter estimation. This function serves as the 'nonlcon' non-linear 
%   constraint function needed by the Optimization Toolbox function FMINCON.
%   The non-linear constraint enforced is the invertibility requirement of
%   moving average polynomial associated with conditional mean specifications
%   with an moving average component.
%
%   [C , CEQ , GC , GCEQ] = garchnlc(Parameters , Series , R , M , P , Q , X)
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
%     in an independent realization (i.e., path). The last row of Series holds 
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
%   C - M-element column vector of non-linear inequality constraints associated
%     with the M roots of the moving average (MA) polynomial. This is an 
%     invertibility constraint that ensures all eigenvalues are inside the 
%     unit circle.
%
%   CEQ - Empty matrix placeholder for future compatibility.
%
%   GC - Empty matrix placeholder for future compatibility.
%
%   GCEQ - Empty matrix placeholder for future compatibility.
%
% Note:
%   This function is needed ONLY for optimization, and is NOT meant to be 
%   called directly. No error checking is performed.
%

% Copyright 1999-2002 The MathWorks, Inc.   
% $Revision: 1.5 $   $ Date: 1998/01/30 13:45:34 $

%
% The non-linear inequality constraints are those associated with 
% the M roots of the moving average (MA) polyomial. The polynomial
% is formed such that the roots computed are the eigenvalues. For 
% an MA noise process to be stable (i.e., invertible) all eigenvalues 
% must lie inside the unit circle of the complex plane. Since the
% 'TolCon' optimization field represents the value by which constraints
% may be violated, the tolerance offset should be set to something 
% larger than TolCon can be expected to be. Setting the tolerance to
% 0.001 implies that eigenvalues can be 1 - 0.001 = 0.999, a highly
% reasonable value.
%

eigenValues =  roots([1 ; Parameters(R+2:R+M+1)]);
tolerance   =  1e-120;
c           =  (abs(eigenValues).^2) - (1 - tolerance);

%
% Return empty matrices as placeholder for future compatibility.
%

ceq  = [];
gc   = [];
gceq = [];
