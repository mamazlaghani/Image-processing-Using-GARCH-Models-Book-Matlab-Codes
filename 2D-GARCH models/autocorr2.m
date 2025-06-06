function [ACF , Lags , bounds] = autocorr2(Series , nLags , nSTDs);
[rows , columns]  =  size(Series);

n           =  length(Series)  % Sample size.
defaultLags =  20;              % BJR recommend about 20 lags for ACFs.

%
% Ensure the number of lags, nLags, is a positive 
% integer scalar and set default if necessary.
%

% Convolution, polynomial multiplication, and FIR digital filtering are
% all the same operation. The FILTER command could be used to compute 
% the ACF (by computing the correlation by convolving the de-meaned 
% Series with a flipped version of itself), but FFT-based computation 
% is significantly faster for large data sets.
%
% The ACF computation is based on Box, Jenkins, Reinsel, pages 30-34, 188.
%

nFFT =  128;
F    =  fft2(Series-mean2(Series),nFFT,nFFT);
F    =  F .* conj(F);
ACF  =  ifft2(F);
ACF  =  ACF((1:(nLags + 1)),(1:(nLags + 1)));         % Retain non-negative lags.
ACF  =  ACF ./ ACF(1);     % Normalize.
ACF  =  real(ACF);

%
% Compute approximate confidence bounds using the Box-Jenkins-Reinsel 
% approach, equations 2.1.13 and 6.2.2, on pages 33 and 188, respectively.
%
sigmaQ  =  sqrt(1/n^2);  
nSTDs=2;
bounds  =  sigmaQ * [nSTDs ; -nSTDs];
Lags    =  [0:nLags]';

%if nargout == 0                     % Make plot if requested.

%
%  Plot the sample ACF.
%
%   lineHandles  =  stem(Lags , ACF , 'filled' , 'r-o');
 %  set   (lineHandles(1) , 'MarkerSize' , 4)
 %  grid  ('on')
 %  xlabel('Lag')
  % ylabel('Sample Autocorrelation')
   %title ('Sample Autocorrelation Function (ACF)')
   %hold  ('on')
%
%  Plot the confidence bounds under the hypothesis that the underlying 
%  Series is really an MA(Q) process. Bartlett's approximation gives
%  an indication of whether the ACF is effectively zero beyond lag Q. 
%  For this reason, the confidence bounds (horizontal lines) appear 
%  over the ACF ONLY for lags GREATER than Q (i.e., Q+1, Q+2, ... nLags).
%  In other words, the confidence bounds enclose ONLY those lags for 
%  which the null hypothesis is assumed to hold. 
%

   %plot([Q+0.5 Q+0.5 ; nLags nLags] , [bounds([1 1]) bounds([2 2])] , '-b');

   
   %plot([0 nLags] , [0 0] , '-k');
   %hold('off')
   %a  =  axis;
   %axis([a(1:3) 1]);

   %clear  ACF  Lags  bounds

   %else

%
%  Re-format outputs for compatibility with the SERIES input. When SERIES is
%  input as a row vector, then pass the outputs as a row vectors; when SERIES
%  is a column vector, then pass the outputs as a column vectors.
%
  % if rowSeries       ACF     =  ACF.';       Lags    =  Lags.';      bounds  =  bounds.';   endend
