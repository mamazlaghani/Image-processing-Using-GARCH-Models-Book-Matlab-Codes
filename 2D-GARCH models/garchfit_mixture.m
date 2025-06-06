function [coefficients,hh, hh1, hh2, e] = garchfit_mixture(spec , y , X)
%GARCHFIT Univariate GARCH process parameter estimation.
%   Given an observed univariate return series, estimate the parameters of a
%   conditional mean specification of ARMAX form and conditional variance 
%   specification of GARCH form. The estimation process infers the innovations 
%   from the return series and fits the model specification to the return 
%   series by maximum likelihood.
%
%   [Coeff,Errors,LLF,Innovations,Sigma,Summary] = garchfit(Series)
%
%   [Coeff,Errors,LLF,Innovations,Sigma,Summary] = garchfit(Spec, Series)
%   [Coeff,Errors,LLF,Innovations,Sigma,Summary] = garchfit(Spec, Series, X)
%
%   garchfit(...)
%
%   Optional Input: Spec, X
%
%   The first calling syntax is strictly a convenience form, modeling a return 
%   series as a constant plus GARCH(1,1) conditionally Gaussian innovations. 
%   For any models beyond this default (yet common) model, a specification 
%   structure, Spec, must be provided.
%
%   The second and third calling syntaxes allow the specification of much more 
%   elaborate models for the conditional mean and variance processes.
%
%   The last calling syntax (with no output arguments) will perform identical
%   estimation procedures as the first three, but will print the iterative 
%   optimization information to the MATLAB command window along with the final
%   parameter estimates and standard errors. It will also produce a tiered plot 
%   of the original return series as well as the innovations (i.e., residuals)
%   inferred, and the corresponding conditional standard deviations.
%
% Inputs:
%   Series - Vector of observations of the underlying univariate return series 
%     of interest. Series is the response variable representing the time series 
%     fit to conditional mean and variance specifications. The last element of 
%     Series holds the most recent observation.
%
% Optional Inputs:
%   Spec - Structure specification for the conditional mean and variance models,
%     and optimization parameters. Spec is a structure with fields created by 
%     calling the function GARCHSET. Type "help garchset" for details.
%
%   X - Time series regression matrix of explanatory variable(s). Typically, X 
%     is a regression matrix of asset returns (e.g., the return series of an 
%     equity index). Each column of X is an individual time series used as an 
%     explanatory variable in the regression component of the conditional mean. 
%     In each column of X, the first row contains the oldest observation and 
%     the last row the most recent. If X is specified, the most recent number 
%     of valid (non-NaN) observations in each column of X must equal or exceed 
%     the most recent number of valid observations in Series. When the number 
%     of valid observations in each column of X exceeds that of Series, only 
%     the most recent observations of X are used. If empty or missing, the 
%     conditional mean will have no regression component.
%
% Outputs:
%   Coeff - Structure containing the estimated coefficients. Coeff is of the 
%     same form as the Spec input structure, which allows other GARCH Toolbox 
%     functions (e.g., GARCHSET, GARCHGET, GARCHSIM, and GARCHPRED) to accept 
%     either Spec or Coeff seamlessly.
%
%   Errors - Structure containing the estimation errors (i.e., the standard 
%     errors) of the coefficients. The fields of Errors are a subset of those 
%     found in Coeff or Spec, and correspond to the coefficient fields only.
%
%   LLF - Optimized log-likelihood objective function value associated with the
%     parameter estimates found in Coeff. Optimization is performed by the 
%     FMINCON function of the MATLAB Optimization Toolbox.
%
%   Innovations - Innovations vector inferred from the input Series. The size
%     of Innovations is the same as the size of Series.
%
%   Sigma - Conditional standard deviation vector corresponding to Innovations.
%     The size of Sigma is the same as the size of Series.
%
%   Summary - Structure of summary information about the optimization process,
%     including convergence information, iterations, objective function calls,
%     active constraints, and the covariance matrix of coefficient estimates.
%
% See also GARCHSET, GARCHSIM, GARCHPRED, GARCHLLFN, FMINCON.

% Copyright 1999-2002 The MathWorks, Inc.   
% $Revision: 1.9 $   $Date: 2002/03/11 19:37:16 $

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
% Check input arguments. The single input case must specify a 
% valid numeric vector of returns. For the multiple input case, 
% the first input must be a specification structure.
%

switch nargin
   case 1
     if isnumeric(spec)
        y    =  spec;
        spec =  garchset;       % Allow a convenience/default model form.
     else
        error(' Observed return series ''Series'' must be specified.');
     end

   case {2 , 3}
     if ~isstruct(spec)
        error(' ''Spec'' must be a structure.');
     end

   otherwise
     error(' Too many inputs specified.');
end

%
% Scrub the observed return series vector y(t).
%

rowY  =  logical(0);     

if prod(size(y)) == length(y)   % Check for a vector (single return series).
   rowY  =  size(y,1) == 1;     % Flag a row vector for outputs.
   y     =  y(:);               % Convert to a column vector.
else
   error(' Observed return series ''Series'' must be a vector.');
end

%
%  The following code segment assumes that missing observations are indicated
%  by the presence of NaN's. Any initial rows with NaN's are removed, and 
%  processing proceeds with the remaining block of contiguous non-NaN rows. 
%  Put another way, NaN's are allowed, but they MUST appear as a contiguous 
%  sequence in the initial rows of the y(t) vector.
%

i1  =  find(isnan(y));
i2  =  find(isnan(diff([y ; zeros(1,size(y,2))]) .* y));

if (length(i1) ~= length(i2)) | any(i1 - i2)
   error(' Only initial observations in ''Series'' may be missing (NaN''s).')
end

if any(sum(isnan(y)) == size(y,1))
   error(' A realization of ''Series'' is completely missing (all NaN''s).')
end

firstValidRow  =  max(sum(isnan(y))) + 1;
y              =  y(firstValidRow:end , :);

%
% Scrub the regression matrix and ensure the observed return series vector y(t) 
% and the regression matrix X(t) have the same number of valid (i.e., non-NaN)
% rows (i.e., impose time index compatibility). During estimation, the innovations
% process e(t) must be inferred from the conditional mean specification, which 
% may include a regression component if desired. In contrast to simulation, 
% estimation of the innovations process e(t) is NOT independent of X. 
%

if (nargin >= 3) & ~isempty(X)

   if prod(size(X)) == length(X)   % Check for a vector.
      X  =  X(:);                  % Convert to a column vector.
   end

%
%  Retain the last contiguous block of non-NaN (i.e, non-missing valued) observations only. 
%
   if any(isnan(X(:)))
      X  =  X((max(find(isnan(sum(X,2)))) + 1):end , :);
   end

   if size(X,1) < size(y,1)
      error(' Regression matrix ''X'' has insufficient number of observations.');
   else
      X  =  X(size(X,1) - (size(y,1) - 1):end , :);    % Retain only the most recent samples.
   end

%
%  Ensure number of regression coefficients (if specified) match number of regressors.
%
   regress  =  garchget(spec , 'Regress');             % Regression coefficients.

   if ~isempty(regress)
      if size(X,2) ~= length(garchget(spec , 'Regress'))
         error(' Number of ''Regress'' coefficients unequal to number of regressors in ''X''.');
      end
   end

else

   X        =  [];   % Ensure X exists.
   regress  =  [];

end

%
% Extract model orders & coefficients.
%

R      =  garchget(spec , 'R');       % Conditional mean AR order.
M      =  garchget(spec , 'M');       % Conditional mean MA order.
P      =  garchget(spec , 'P');       % Conditional variance order for lagged variances.
Q      =  garchget(spec , 'Q');       % Conditional variance order for lagged squared residuals.

C      =  garchget(spec , 'C');       % Conditional mean constant.
AR     =  garchget(spec , 'AR');      % Conditional mean AR coefficients.
MA     =  garchget(spec , 'MA');      % Conditional mean MA coefficients.

K      =  garchget(spec , 'K');       % Conditional variance constant.
GARCH  =  garchget(spec , 'GARCH');   % Conditional variance coefficients for lagged variances.
ARCH   =  garchget(spec , 'ARCH');    % Conditional variance coefficients for lagged squared residuals.

%
% Set a 'Display' flag to determine if any warnings or 
% messages should be printed to the screen to alert users.
%

D            =  garchget(spec , 'Display');
DisplayFlag  =  strcmp(D(~isspace(D)) , 'on');

%
% Generate initial parameter estimates if necessary. Note that the code below
% adopts an 'all or nothing' approach for model specifications. That is, if a
% user wants to provide parameter values (either as (1) Initial guesses for 
% refinement during the optimization process, or as (2) Equality constraints
% in which some of the parameters are held fixed and some are refined), the 
% user MUST provide a complete specification. The only flexibility regarding
% model specification is that the conditional mean and variance models are 
% treated separately. For example, a user may provide initial values for the 
% conditional mean parameters, and let the code provide initial values for
% conditional variance parameters.
%

%
% Initialize conditional mean parameter estimates if necessary. If an
% incomplete conditional mean specification is found (i.e., the user has
% NOT explicitly set ALL required coefficients for a given conditional mean 
% model), then compute initial estimates and over-write any pre-existing 
% parameters in the incomplete specification. This is an 'all or nothing' 
% approach for model specification.
%

Flag1                  =  logical(1);  % Initialize to a complete mean specification.
unconditionalVariance  =  [];          % Make sure it exists.

if isempty(X)

   if isempty(C) | (isempty(AR) & (R > 0)) | (isempty(MA) & (M > 0))
%
%     General ARMA conditional mean with no regression component.
%
      [AR , MA , C , unconditionalVariance]  =  arma0(y , R , M);

      Flag1  =  logical(0);   % Indicate an INCOMPLETE mean specification.

   end

else

   if M == 0        % Check for MA terms.

      if isempty(C) | (isempty(AR) & (R > 0)) | isempty(regress)
%
%        General ARX conditional mean model with no MA terms. Initial
%        estimates can be generated by a simple OLS regression.
%
         yLag               =  lagmatrix(y , [1:R]);
         yLag(isnan(yLag))  =  mean(y);

         [QQ , RR]  =  qr([ones(size(X,1),1)  yLag  X] , 0);
         b          =  RR \ (QQ'*y);
         residuals  =  y - [ones(size(X,1),1)  yLag  X]*b;

         C       =  b(1);
         AR      =  b(2:R+1);
         MA      =  [];
         regress =  b(R+2:end);

         unconditionalVariance  =  var(residuals,1);
%
%        Ensure stationarity of AR process.
%
         if any(abs(roots([1 ; -AR(:)])) >= 1)
            AR(:)  =  0;
         end

         Flag1  =  logical(0);            % Indicate an INCOMPLETE mean specification.

      end

   else

      if isempty(C) | (isempty(AR) & (R > 0)) | isempty(MA) | isempty(regress)
%
%        General ARMAX conditional mean model. 
%
         C       =  0; 
         AR      =  zeros(R,1);
         MA      =  zeros(M,1); 
         regress =  zeros(size(X,2),1);

         unconditionalVariance  =  1e-4;  % Default assumption for daily returns.

         Flag1  =  logical(0);            % Indicate an INCOMPLETE mean specification.

      end

   end

end

%
% Create initial conditional variance parameter estimates if necessary.
% As for the conditional mean above, this is an 'all or nothing' approach
% for specification of the conditional variance.
%

Flag2  =  logical(1);           % Initialize flag to indicate a complete variance specification.

if isempty(K) | (isempty(GARCH) & (P > 0)) | (isempty(ARCH) & (Q > 0))

   [K , GARCH , ARCH] = garch0(P , Q , unconditionalVariance);

   Flag2  =  logical(0);        % Indicate an INCOMPLETE variance specification.

end

%
% Extract equality constraint information for coefficients.
%

Fix  =  zeros(1 + R + M + size(X,2) + 1 + P + Q , 1);

if Flag1
%
%  A complete conditional mean specification was provided.
%
   FixC       =  garchget(spec , 'FixC');
   FixAR      =  garchget(spec , 'FixAR');
   FixMA      =  garchget(spec , 'FixMA');
   FixRegress =  garchget(spec , 'FixRegress');

   if isempty(FixC)       , FixC       = 0;                  end
   if isempty(FixAR)      , FixAR      = zeros(R,1);         end
   if isempty(FixMA)      , FixMA      = zeros(M,1);         end
   if isempty(FixRegress) , FixRegress = zeros(size(X,2),1); end

   Fix(1:(1 + R + M + size(X,2))) = [FixC ; FixAR(:) ; FixMA(:) ; FixRegress(:)];

end

if Flag2
%
%  A complete conditional variance specification was provided.
%
   FixK      =  garchget(spec , 'FixK');
   FixGARCH  =  garchget(spec , 'FixGARCH');
   FixARCH   =  garchget(spec , 'FixARCH');

   if isempty(FixK)       , FixK      = 0;                  end
   if isempty(FixGARCH)   , FixGARCH  = zeros(P,1);         end
   if isempty(FixARCH)    , FixARCH   = zeros(Q,1);         end

   Fix((2 + R + M + size(X,2)):end) = [FixK ; FixGARCH(:) ; FixARCH(:)];

end

%
% NOTES:
% Let y(t) = return series of interest (assumed stationary)
%     e(t) = innovations of the model noise process (assumed invertible)
%     h(t) = conditional variance of the innovations process e(t)
%
% The initial parameter guess vector 'x0' input to FMINCON and the estimated 
% parameter vector 'coefficient' output from FMINCON are formatted exactly as
% the parameters would be read from the recursive difference equations when 
% solving for the current values of the y(t) and h(t) time series.
%
% For example, consider the following equations for the conditional mean 
% and variance of an ARMAX(R=2,M=2,nX=1)/GARCH(P=2,Q=2) composite model:
%
%   y(t) =  1.3 + 0.5y(t-1) - 0.8y(t-2) + e(t) - 0.6e(t-1) + 0.08e(t-2) + 1.2X(t)
%   h(t) =  0.5 + 0.2h(t-1) + 0.1h(t-2) + 0.3e(t-1)^2  +  0.4e(t-2)^2
%
% The above equations are examples of the following general 
% ARMAX(R,M,nX)/GARCH(P,Q) form:
%
%   y(t) =  C + AR(1)y(t-1) + ... + AR(R)y(t-R) + e(t) 
%             + MA(1)e(t-1) + ... + MA(M)e(t-M) + B(1)X(t,1) + ... + B(nX)X(t,nX)
%
%   h(t) =  K + GARCH(1)h(t-1)   + ... + GARCH(P)h(t-P) 
%             +  ARCH(1)e(t-1)^2 + ... +  ARCH(Q)e(t-Q)^2
%
% For the example listed above, the coefficient vector would be formatted as
%
%   Coefficient = [ C     AR(1:R)     MA(1:M)  B(1:nX)  K   GARCH(1:P)   ARCH(1:Q)]'
%               = [1.3   0.5 -0.8   -0.6 0.08    1.2   0.5    0.2 0.1     0.3 0.4 ]'
%
% Notice that the coefficient of e(t) in the conditional mean equation is
% defined to be 1, and is NOT included in the vector because it is not estimated.
%

%
% Get the probability distribution of the innovations process e(t) 
% and call the appropriate log-likelihood objective function.
%





      x0  =  [C ; AR(:) ; MA(:) ; regress(:) ; K ; GARCH(:) ; ARCH(:) ; K ; GARCH(:) ; ARCH(:); .9];   % Initial guess.
size(x0);
      %

%     Set lower bounds constraints. 
%
%     The parameters of the conditional mean, in theory, have no lower bounds. 
%     However, for robustness, set the mean constant (C) and the ARMA(R,M) 
%     parameters to -100 (i.e., an unrealistically large, yet totally finite 
%     number). The coefficients of the regression component of the mean are 
%     set to -infinity to reflect the fact that the origin and generating 
%     mechanism of a regression is completely unknown. This is consistent 
%     with the 'hands-off' approach regarding regression.
%
%     The conditional variance constant (K) must be positive, so it's set to 
%     a very small (yet positive) value; the conditional variance coefficients 
%     (GARCH(i), ARCH(i)) are constrained to be non-negative, so they're set 
%     to zero.
%
      infinity     =  inf;
      hundred      =  100;
      lowerBounds  =  [-hundred(ones(1+R+M,1)) ; -infinity(ones(size(X,2),1)) ; 1e-10 ; zeros(P+Q,1); 1e-10 ; zeros(P+Q,1);0];
      
      upperBounds  =  [infinity(ones(1+R+M+1+P+Q+1+P+Q,1)); 1];

      size(upperBounds);
%
%     Set linear inequality of the covariance-stationarity constraint of 
%     the conditional variance, Ax <= b. Since the conditional variance 
%     (GARCH(i), ARCH(i)) parameters are constrained to be non-negative, 
%     the covariance-stationarity constraint is just a summation constraint.
%     Also, adjust the summation constraint, b, to reflect a tolerance 
%     offset from a fully integrated conditional variance condition (i.e.,
%     an IGARCH process).
%

      if ((P + Q) > 0)
         A  =  [zeros(1,1+R+M+size(X,2)+2+P+Q)  ones(1,P)  ones(1,Q) 0];
         size(A);
         b  = 0.99999800000000;
      else
         A  =  [];
         b  =  [];
      end
%
%     Set any linear equality constraints.
%
      Aeq=[];
      beq=[];
%
%     Perform constrained non-linear optimization. 
%

     [coefficients , LLF, exitFlag , output , lambda] =  fmincon('garchllfn_mixture'  , x0 , A  , b, Aeq , beq ,lowerBounds ,upperBounds ,'garchnlc_mixture', [] ,y , R , M, P , Q,  X);
%
%     Negate objective function value to compensate for FMINCON.
%
      LLF  =  -LLF;
%
%     Over-write all GARCH constraint-violating parameters that are less than 
%     zero. This will, occasionally, occur because FMINCON may violate constraints 
%     ever so slightly. Also, the constant term (K) of the conditional variance 
%     equation must be positive, so enforce this if necessary.
%
      [LLF , G , H , e , hh, hh1, hh2] = garchllfn_mixture(coefficients , y , R, M , P ,Q , X);
