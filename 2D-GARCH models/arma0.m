function  [AR, MA, constant, variance] = arma0(y , R , M);
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
   correlation =  autocorr(y , R + M);      % auto-correlation sequence.
   variance    =  var(y,1);
   covariance  =  correlation * variance ;   % auto-covariance sequence.

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
      i(i <= 0)  =  i(i <= 0) + 2;     % covariance(k) = covariance(-k)

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

