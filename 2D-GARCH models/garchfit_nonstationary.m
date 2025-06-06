function [coefficients, errors, LLF, innovations, sigma, summary] = garchfit_nonstationary(spec , y , X)
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

distribution  =  garchget(spec , 'Distribution');
distribution  =  distribution(~isspace(distribution));

switch upper(distribution)

   case 'GAUSSIAN'

      x0  =  [C ; AR(:) ; MA(:) ; regress(:) ; K ; GARCH(:) ; ARCH(:)];   % Initial guess.
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
      lowerBounds  =  [-hundred(ones(1+R+M,1)) ; -infinity(ones(size(X,2),1)) ; 1e-10 ; zeros(P+Q,1)];
      upperBounds  =  [];

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
         %A  =  [zeros(1,1+R+M+size(X,2)+1)  ones(1,P)  ones(1,Q)];
         %b  =  1  -  2*optimget(spec.Optimization , 'TolCon', 1e-6);
         A=[];
         b=[];
      else
         A  =  [];
         b  =  [];
      end
%
%     Set any linear equality constraints.
%
      if any(Fix)
         i    =  find(Fix); 
         Aeq  =  [zeros(length(i) , 1 + R + M + size(X,2) + 1 + P + Q)];

         for j = 1:length(i)
             Aeq(j,i(j))  =  1;
         end

         beq  =  x0(logical(Fix));

      else
         Aeq  =  [];
         beq  =  [];
      end
%
%     Perform constrained non-linear optimization. 
%

      [coefficients , LLF    , ...
       exitFlag     , output , lambda] =  fmincon('garchllfn'  , x0 , A  , b       , Aeq , beq , ...
                                                   lowerBounds , upperBounds       , ...
                                                  'garchnlc'   , spec.Optimization , ...
                                                   y , R , M , P , Q , X);
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
      varianceCoefficients =  coefficients((2 + R + M + size(X,2)):end);

      varianceCoefficients(varianceCoefficients <  0) =  0;              % All variance-related coefficients.
      varianceCoefficients(1) = max(varianceCoefficients(1) , realmin);  % Constant term 'K'.

      coefficients  =  [coefficients(1:(1 + R + M + size(X,2))) ; varianceCoefficients(:)];

      MLEparameters =  coefficients(:);  % Save MLE parameter values in vector form.
%
%     Extract the parameter estimates & equality constraints.
%
      C          =  coefficients(1);
      AR         =  coefficients(2:R+1);
      MA         =  coefficients(R+2:R+M+1);
      regress    =  coefficients(R+M+2:R+M+size(X,2)+1);

      K          =  coefficients(R+M+size(X,2)+2:R+M+size(X,2)+2);
      GARCH      =  coefficients(R+M+size(X,2)+3:R+M+size(X,2)+2+P);
      ARCH       =  coefficients(R+M+size(X,2)+3+P:R+M+size(X,2)+2+P+Q);

      FixC       =  Fix(1);
      FixAR      =  Fix(2:R+1);
      FixMA      =  Fix(R+2:R+M+1);
      FixRegress =  Fix(R+M+2:R+M+size(X,2)+1);

      FixK       =  Fix(R+M+size(X,2)+2:R+M+size(X,2)+2);
      FixGARCH   =  Fix(R+M+size(X,2)+3:R+M+size(X,2)+2+P);
      FixARCH    =  Fix(R+M+size(X,2)+3+P:R+M+size(X,2)+2+P+Q);

%
%     Test for stationarity & invertibility of the ARMA model (if any).
%
      summary.warning = 'No Warnings';

      if any(abs(roots([1 ; -AR(:)])) >= 1)
         summary.warning = 'ARMA Model is Not Stationary/Invertible';
         if DisplayFlag
            warning('Auto-Regressive Polynomial is Non-Stationary.')
         end
      end

      if any(abs(roots([1 ; MA(:)])) >= 1)
         summary.warning = 'ARMA Model is Not Stationary/Invertible';
         if DisplayFlag
            warning('Moving-Average Polynomial is Non-Invertible.')
         end
      end
%
%     Now pack the data into the output COEFFICIENTS structure. Note that 
%     the COEFFICIENTS output structure is of the same form as the SPEC 
%     input structure. This allows GARCHSET, GARCHGET, GARCHSIM, and 
%     GARCHPRED to accept either SPEC or COEFFICIENTS seamlessly. 
%
%     Strictly speaking, GARCHSET should be used to make the following
%     assignment, but will error-out if the stationarity/invertibility
%     constraints are violated. Simple assignment allows for graceful 
%     termination with a warning.
%
      coefficients            =  spec;
      coefficients.C          =  C;
      coefficients.AR         =  AR(:)';
      coefficients.MA         =  MA(:)';
      coefficients.Regress    =  regress(:)';
      coefficients.K          =  K;
      coefficients.GARCH      =  GARCH(:)';
      coefficients.ARCH       =  ARCH(:)';

      coefficients.FixC       =  FixC;
      coefficients.FixAR      =  FixAR(:)';
      coefficients.FixMA      =  FixMA(:)';
      coefficients.FixRegress =  FixRegress(:)';
      coefficients.FixK       =  FixK;
      coefficients.FixGARCH   =  FixGARCH(:)';
      coefficients.FixARCH    =  FixARCH(:)';
%
%     Update the 'comment' field to reflect a regression component (only if auto-generated).
%
      comment =  garchget(coefficients , 'comment');

      if length(findstr(comment, char(0))) == 2
         pOpen  =  findstr(comment, '(');
         pClose =  findstr(comment, ')');
         if ~isempty(pOpen) & ~isempty(pClose) 
            commas  =  findstr(comment(pOpen(1):pClose(1)) , ',');
            if length(commas) == 1
               coefficients.Comment =  [comment(1:pClose(1)-1) ',' num2str(size(X,2)) comment(pClose(1):end)];
            elseif length(commas) == 2
               coefficients.Comment =  [comment(1:(pOpen(1) + commas(2)-1)) num2str(size(X,2)) comment(pClose(1):end)];
            end
        end
      end

      if length(coefficients.AR)      == 0 , coefficients.AR          =  []; end   % Just for aesthetics.
      if length(coefficients.MA)      == 0 , coefficients.MA          =  []; end
      if length(coefficients.Regress) == 0 , coefficients.Regress     =  []; end
      if length(coefficients.GARCH)   == 0 , coefficients.GARCH       =  []; end
      if length(coefficients.ARCH)    == 0 , coefficients.ARCH        =  []; end

      if sum(coefficients.FixC)       == 0 , coefficients.FixC        =  []; end   % Just for aesthetics.
      if sum(coefficients.FixAR)      == 0 , coefficients.FixAR       =  []; end
      if sum(coefficients.FixMA)      == 0 , coefficients.FixMA       =  []; end
      if sum(coefficients.FixRegress) == 0 , coefficients.FixRegress  =  []; end
      if sum(coefficients.FixK)       == 0 , coefficients.FixK        =  []; end
      if sum(coefficients.FixGARCH)   == 0 , coefficients.FixGARCH    =  []; end
      if sum(coefficients.FixARCH)    == 0 , coefficients.FixARCH     =  []; end

%
%     Compute the variance-covariance matrix of the parameter 
%     estimates and extract the standard errors of the estimation.
%
      if (nargout >= 2) | (nargout == 0)

         covarianceMatrix =  varcov(MLEparameters , y , R , M , P , Q , X , Fix);
         standardErrors   =  sqrt(diag(covarianceMatrix))';

         errors.C         =  standardErrors(1);
         errors.AR        =  standardErrors(2:R+1);
         errors.MA        =  standardErrors(R+2:R+M+1);
         errors.Regress   =  standardErrors(R+M+2:R+M+size(X,2)+1);

         errors.K         =  standardErrors(R+M+size(X,2)+2:R+M+size(X,2)+2);
         errors.GARCH     =  standardErrors(R+M+size(X,2)+3:R+M+size(X,2)+2+P);
         errors.ARCH      =  standardErrors(R+M+size(X,2)+3+P:R+M+size(X,2)+2+P+Q);

         if length(errors.AR)      == 0 , errors.AR       =  []; end   % Just for aesthetics.
         if length(errors.MA)      == 0 , errors.MA       =  []; end
         if length(errors.Regress) == 0 , errors.Regress  =  []; end
         if length(errors.GARCH)   == 0 , errors.GARCH    =  []; end
         if length(errors.ARCH)    == 0 , errors.ARCH     =  []; end

      end
%
%     Return the innovations and conditional standard deviation vectors if requested.
%
      if (nargout >= 4) | (nargout == 0)
         [innovations , sigma]  =  garchinfer(coefficients , y , X);
      end
%
%     Return summary information if requested.
%
      if (nargout >= 6) | (nargout == 0)

         if exitFlag == 0
            summary.converge  =  'Maximum Function Evaluations or Iterations Reached';
         elseif exitFlag < 0
            summary.converge  =  'Function Did NOT Converge';
         elseif exitFlag > 0
            summary.converge  =  'Function Converged to a Solution';
         end
 
         summary.covMatrix      =  covarianceMatrix;
         summary.iterations     =  output.iterations;
         summary.functionCalls  =  output.funcCount;
%
%        Flag any boundary constraints enforced EXCEPT LINEAR EQUALITY CONSTRAINTS
%        specifically requested by the user. Whenever any constraints (excluding
%        linear equalities) are imposed, the log-likelihood function will probably NOT
%        be approximately quadratic at the solution. In this case, the standard errors
%        of the parameter estimates are unlikely to be accurate. However, if the ONLY
%        constraints are the linear equalities imposed by the user, then the resulting
%        log-likelihood value LLF may still be useful for post-fit assessment and
%        inference tests, such as likelihood ratio tests (LRT's).
%

         TolCon = optimget(spec.Optimization , 'TolCon', 1e-6);

         if (norm([lambda.lower(:) ; lambda.upper(:) ; lambda.ineqlin(:) ; lambda.ineqnonlin(:)] , 1) > TolCon)
            summary.constraints  =  'Boundary Constraints Active; Errors may be Inaccurate';
            if DisplayFlag
               warning('Boundary Constraints Active; Standard Errors may be Inaccurate.')
            end
         else
            summary.constraints  =  'No Boundary Constraints';
         end

      end

   otherwise

      error(' Distribution of innovations must be ''Gaussian''.')

end

%
% Re-format outputs for compatibility with the SERIES input. When 
% SERIES is input as a single row vector, then pass the outputs 
% as a row vectors. 
%

if rowY & (nargout >= 4)
   innovations  =  innovations(:).';
   sigma        =  sigma(:).';
end

%
% Perform the default no-output action: 
%
%  (1) Print the parameter estimates to the screen, and 
%  (2) Display the estimated residuals, conditional standard
%      deviations, and input raw return series.
%

if nargout == 0

   garchdisp(coefficients , errors);
   garchplot(innovations  , sigma , y);

   disp(' ')
   fprintf('  Log Likelihood Value: %f\n\n' , LLF)

   clear coefficients     % Suppress unexpected printing to the command window.

end

%
%   * * * * *  Helper function for initial GARCH guesses.  * * * * *
%

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


%
%   * * * * *  Helper function for initial ARMA guesses.  * * * * *
%

function [AR, MA, constant, variance] = arma0(y , R , M)
%ARMA0 Initial parameter estimates of univariate ARMA processes.
%   Compute initial estimates of the auto-regressive (AR) and moving average 
%   (MA) coefficients of a stationary/invertible univariate ARMA time series. 
%   Estimates of the model constant and the variance of the innovations noise 
%   process are also provided. The purpose of this function is to provide 
%   initial coefficient estimates suitable for further refinement via maximum 
%   likelihood estimation.
%
%   [AR , MA , Constant , Variance] = arma0(Series)
%   [AR , MA , Constant , Variance] = arma0(Series , R , M)
%
%   Optional Inputs: R, M
%
% Input:
%   Series - Return series vector of interest. Series is a dependent stochastic 
%     process assumed to follow a model specification of general ARMA(R,M) form. 
%     The first element contains the oldest observation and the last element 
%     the most recent. 
%
% Optional Inputs:
%   R - Auto-regressive order of an ARMA(R,M) model. R is a non-negative integer
%     scalar. If empty or missing, R = 0.
%
%   M - Moving average order of an ARMA(R,M) model. M is a non-negative integer 
%     scalar. If empty or missing, M = 0.
%
% Outputs:
%   AR - R by 1 vector of initial estimates of the auto-regressive parameters 
%     of an ARMA(R,M) model. The first element of AR is an estimate of the 
%     coefficient of the first lag of Series, the second element is an estimate
%     of the coefficient of the second lag of Series, and so forth.
%      
%   MA - M by 1 vector of initial estimates of the moving average parameters of 
%     an ARMA(R,M) model. The first element of MA is an estimate of the 
%     coefficient of the first lag of the innovations noise process, the second
%     element is an estimate of the coefficient of the second lag of the 
%     innovations noise process, and so forth.
%
%   Constant - Estimate of the constant included in a general ARMA model. 
%
%   Variance - Estimate of the unconditional variance of the innovations noise
%     process. Equivalently, it may also be viewed as the variance estimate of
%     the white noise innovations under the assumption of homoskedasticity.
%
% See also GARCHSET, GARCHGET, GARCHSIM, GARCHFIT.

%
% References:
%   Box, G.E.P., Jenkins, G.M., Reinsel, G.C., "Time Series Analysis: 
%     Forecasting and Control", 3rd edition, Prentice Hall, 1994.
%   Hamilton, J.D., "Time Series Analysis", Princeton University Press, 1994.
%

%
% Check & scrub the observed return series matrix y(t).
%

if (nargin < 1) | isempty(y)
   error(' Observed return series ''Series'' must be specified.');
else
   if prod(size(y)) == length(y)   % Check for a vector (single return series).
      y  =  y(:);                  % Convert to a column vector.
   end
end

if (nargin >= 2) & ~isempty(R)
   if length(R) ~= 1
      error(' AR order ''R'' must be a scalar.');
   end
   if any(round(R) - R) | any(R < 0)
      error(' AR order ''R'' must be a non-negative integer.')
   end
else
   R  =  0;
end

if (nargin >= 3) & ~isempty(M)
   if length(M) ~= 1
      error(' MA order ''M'' must be a scalar.');
   end
   if any(round(M) - M) | any(M < 0)
      error(' MA order ''M'' must be a non-negative integer.')
   end
else
   M  =  0;
end

%
% Check for the simple case of no ARMA model. In this 
% case just compute the sample mean and variance and exit.
%

if (R + M) == 0
   AR        =  [];
   MA        =  [];
   constant  =  mean(y);
   variance  =  var(y,1);
   return
end

%
% Estimate the AR coefficients of a general ARMA(R,M) model.
%

if R > 0

%
%  Compute the auto-covariance sequence of the y(t) process. Note that the
%  variance (i.e, zeroth-lag auto-covariance) is found in the first element.
%
   correlation =  autocorr(y , R + M);         % auto-correlation sequence.
   variance    =  var(y,1);
   covariance  =  correlation * variance;      % auto-covariance sequence.

%
%  In each case below, the matrix 'C' of covariances is 
%  that outlined in equation A6.2.1 (page 220) of BJR.
%
   if M > 0
%
%     For ARMA processes, the matrix C of covariances derived from the 
%     estimated auto-covariance sequence is Toeplitz, but non-symmetric.
%     The AR coefficients are then found by solving the modified Yule-Walker
%     equations.
%
      i          =  [M+1:-1:M-R+2]; 
      i(i <= 0)  =  i(i <= 0) + 2;       % covariance(k) = covariance(-k)

      C  =  toeplitz(covariance(M+1:M+R) , covariance(i));

      if R == 1
         AR  =  covariance(M+2:M+R+1) / C;
      else
         AR  =  C \ covariance(M+2:M+R+1);
      end

   else

      if R == 1
         AR  =  correlation(2);
      else
%
%        For AR processes, the matrix C of covariances derived from the 
%        estimated auto-covariance sequence is Toeplitz and symmetric.
%        The AR coefficients are found by solving the Yule-Walker equations.
%
         C   =  toeplitz(covariance(1:R));
         AR  =  C \ covariance(2:R+1);

      end

   end

%
%  Ensure the AR process is stationary. If it's not stationary, then set all
%  ARMA coefficients to 0. This ensures the subsequent optimization will
%  begin with a stationary/invertible ARMA model for the conditional mean.
%

   eigenValues =  roots([1 ; -AR(:)]);

   if any(abs(eigenValues) >= 1)

      AR        =  zeros(R , 1);
      MA        =  zeros(M , 1);
      constant  =  0;
      variance  =  1e-4;   % A decent assumption for daily returns.
      return

   end

else

   AR  =  [];

end

%
% Filter the ARMA(R,M) input series y(t) with the estimated AR coefficients 
% to obtain a pure MA process. If the input moving-average model order M is
% zero (M = 0), then the filtered output is really just a pure innovations 
% process (i.e., an MA(0) process); in this case the innovations variance 
% estimate is just the sample variance of the filtered output. If M > 0, then
% compute the auto-covariance sequence of the MA process and continue.
%

x        =  filter([1 -AR'] , 1 , y);
constant =  mean(x);

if M == 0
   variance =  var(x,1);
   MA       =  [];
   return
end

c  =  autocorr(x , M) * var(x,1);    % Covariance of an MA(M) process.

%
% Estimate the variance of the white noise innovations process e(t)
% and the MA coefficients of a general ARMA(R,M) model. The method of
% computation is that outlined in equation A6.2.4 (page 221) of BJR.
%

MA      =  zeros(M , 1);  % Initialize MA coefficients.
MA1     =  ones (M , 1);  % Saved MA coefficients from previous iteration.
counter =  1;             % Iteration counter.
tol     =  0.05;          % Convergence tolerance.

while ((norm(MA - MA1) > tol) & (counter < 100))

    MA1  =  MA;
%
%   Estimate the variance of the innovations process e(t).
%
    variance  =  c(1) /([1 ; MA]'* [1 ; MA]);

    if abs(variance) < tol 
       break
    end

    for j = M:-1:1
%
%       Now estimate the moving-average coefficients. Note that 
%       the MA coefficients are the negative of those appearing 
%       in equation A6.2.4 (page 221) of BJR. This is due to the
%       convention of entering coefficient values, via GARCHSET,
%       exactly as the equation would be written.
%
        MA(j)  =  [c(j+1) ; -MA(1:M-j)]' * [1/variance ; MA(j+1:M)];

    end

    counter  =  counter  +  1;

end

%
% Test for invertibility of the noise model.
%

eigenValues  =  roots([1 ; MA]);

if any(abs(eigenValues) >= 1) | any(isinf(eigenValues)) | any(isnan(eigenValues))

   MA        =  zeros(M , 1);
   variance  =  1e-4;   % A decent assumption for daily returns.

end


%
%   * * * * *  Helper function for variance-covariance.  * * * * *
%

function covarianceMatrix = varcov(p0 , y , R , M , P , Q , X , Fix)
%VARCOV Covariance matrix of maximum likelihood parameter estimates.
%   This function uses the outer-product method to compute the error 
%   covariance matrix of parameters estimated by maximum likelihood. It 
%   is valid for maximum likelihood estimation only, and assumes that 
%   the peak of the log-likelihood objective function has been found 
%   within the interior of the allowable parameter space (i.e., no 
%   boundary constraints have been actively enforced).
%
%   V = varcov(Parameters , Series , R , M , P , Q , X , Fix)
%
% Inputs:
%   Parameters - Column vector of optimized maximum likelihood parameter 
%     estimates associated with fitting conditional mean and variance 
%     specifications to an observed return series Series. The conditional 
%     mean contributes the first (1 + R + M + nX) parameters (nX is the 
%     number of explanatory variables included in the regression component 
%     of the conditional mean); the conditional variance contributes the 
%     remaining (1 + P + Q) parameters. The length of Parameters is thus 
%     (2 + R + M + nX + P + Q) for a Gaussian log-likelihood function.
%
%   Series - Column vector of observations of the underlying univariate return
%     series of interest. Series is the response variable representing the time
%     series fit to conditional mean and variance specifications. The last row 
%     of Series holds the most recent observation of each realization.
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
%   Fix - Boolean (0,1) vector the same length as Parameters. 0's indicate 
%     that the corresponding Parameter element has been estimated; 1's indicate 
%     that the corresponding Parameter element has been held fixed.
%
% Outputs:
%   V - Variance-covariance matrix of the estimation errors associated with
%     parameter estimates obtained by maximum likelihood estimation (MLE). The
%     standard errors of the estimates are the square root of the diagonal 
%     elements. V is computed by the 'outer-product' method.
%
% Note: 
%   This is an internal helper function for GARCHFIT. No error checking is 
%   performed.
%

% References:
%   Hamilton, J.D., "Time Series Analysis", Princeton University 
%     Press, 1994, pages 142-144, 660-661.
%

delta  =  1e-10;  % Offset for numerical differentiation.

%
% Evaluate the log-likelihood objective function at the final MLE 
% parameter estimates. In contrast to the optimization which is 
% interested in the scalar-valued objective function 'LLF', here 
% we are interested in the individual log-likelihood components 
% for each observation of y(t). The relationship between them
% 
%                      LLF = -sum(g0)
%

[LLF,G,H,residuals,sigma] =  garchllfn(p0 , y , R , M , P , Q , X);
g0                        = -0.5 * ( log(2*pi*(sigma.^2)) + (residuals./sigma).^2 );

%
% Initialize the perturbed parameter vector and the scores matrix. 
% For 'T' observations in y(t) and 'nP' parameters estimated via
% maximum likelihood, the scores array is a T-by-nP matrix.
%

pDelta  =  p0;
scores  =  zeros(length(y) , length(p0));

for j=1:length(p0)

   if ~Fix(j)

      pDelta(j)   =  p0(j) * (1+delta);
      dp          =  delta * p0(j);
%
%     Trap the case of a zero parameter value (i.e., p0(j) = 0).
%
      if dp == 0
         dp        =  delta;
         pDelta(j) =  dp;
      end

      [LLF,G,H,residuals,sigma] = garchllfn(pDelta , y , R , M , P , Q , X);

      gDelta      =  -0.5 * ( log(2*pi*(sigma.^2)) + (residuals./sigma).^2 );

      scores(:,j) =  (g0 - gDelta) / dp;
      pDelta(j)   =  p0(j);

   end

end

%
% Invert the outer product of the scores matrix to get the 
% variance-covariance matrix of the MLE parameter estimates.
%

covarianceMatrix  =  pinv(scores'*scores);
