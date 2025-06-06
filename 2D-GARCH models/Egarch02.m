function [K , GARCH , ARCH] = Egarch02(p1, p2 , q1, q2 , unconditionalVariance)
GARCH=0*ones((p1+1)*(p2+1)-1,1);
ARCH=0*ones((q1+1)*(q2+1)-1,1);
GARCH =  GARCH(:);
ARCH  =  ARCH(:);
if isempty(unconditionalVariance) |  (unconditionalVariance <= 0)
   K  =  1e-5;   % A decent assumption for daily returns.
else
   K  =  unconditionalVariance * (1 - sum([GARCH ; ARCH]));
end
