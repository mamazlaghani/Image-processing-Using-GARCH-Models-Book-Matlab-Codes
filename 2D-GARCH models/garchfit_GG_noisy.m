function [coefficients,h]=garchfit_GG_noisy(spec , y , X);
%
%y1=y(2:end);van=y(1);
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
D            =  garchget(spec , 'Display');
DisplayFlag  =  strcmp(D(~isspace(D)) , 'on');
Flag1                  =  logical(1);  % Initialize to a complete mean specification.
unconditionalVariance  =  [];          % Make sure it exists.

van=y(1);
y=y(2:end);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
x0  =  [C ; AR(:) ; MA(:) ; regress(:) ; K ; GARCH(:) ; ARCH(:);.4,;.6];   % Initial guess.
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

infinity     =  inf;
      hundred      =  100;
      %lowerBounds  =  [-hundred(ones(1+R+M,1)) ; -infinity(ones(size(X,2),1)) ; 1e-10 ; zeros(P+Q+1,1);1];
      lowerBounds  =  [-hundred(ones(1+R+M,1)) ; -infinity(ones(size(X,2),1)) ; 1e-10 ; zeros(P+Q,1);.001;1];
      %upperBounds = []
      upperBounds  =  [infinity(ones(1+R+M,1));infinity(ones(size(X,2),1));infinity(ones(P+Q+2,1));2 ];
      % Set linear inequality of the covariance-stationarity constraint of 
%     the conditional variance, Ax <= b. Since the conditional variance 
%     (GARCH(i), ARCH(i)) parameters are constrained to be non-negative, 
%     the covariance-stationarity constraint is just a summation constraint.
%     Also, adjust the summation constraint, b, to reflect a tolerance 
%     offset from a fully integrated conditional variance condition (i.e.,
%     an IGARCH process).
      if ((P + Q) > 0)
         A  =  [zeros(1,1+R+M+size(X,2)+1)  ones(1,P)  ones(1,Q) zeros(1,2)];
         b  =  1  -  2*optimget(spec.Optimization , 'TolCon', 1e-6);
      else
         A  =  [];
         b  =  [];
      end
 if any(Fix)
         i    =  find(Fix); 
         Aeq  =  [zeros(length(i) , 1 + R + M + size(X,2) + 1 + P + Q+2)];

         for j = 1:length(i)
             Aeq(j,i(j))  =  1;
         end

         beq  =  x0(logical(Fix));

      else
         Aeq  =  [];
         beq  =  [];
      end
  end
%
[coefficients , LLF, exitFlag , output , lambda] =  fmincon('garchllfn_GG_noisy1'  , x0 , A  , b, Aeq , beq ,lowerBounds ,upperBounds ,[] , spec.Optimization,y , R , M , P , Q ,van, X);


[LLF , G , H , e , h] = garchllfn_GG_noisy1(coefficients , y , R , M , P , Q ,van, X);






 %LLF  =  -LLF;
%
%     Over-write all GARCH constraint-violating parameters that are less than 
%     zero. This will, occasionally, occur because FMINCON may violate constraints 
%     ever so slightly. Also, the constant term (K) of the conditional variance 
%     equation must be positive, so enforce this if necessary.
%
  %    varianceCoefficients =  coefficients((2 + R + M + size(X,2)):end);

     