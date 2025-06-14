function [coefficients,hh,e] = Egarchfitt2(spec,y,X);

r1=spec(1);
r2=spec(2);
m1=spec(3);
m2=spec(4);
p1=spec(5);
p2=spec(6);
q1=spec(7);
q2=spec(8);

C=[];
regress  =  [];
AR=[];
MA=[];
GARCH=[];
ARCH=[];
K=[];
unconditionalVariance  =  [];          % Make sure it exists.

if isempty(X)

   if isempty(C) | (isempty(AR) & (r1 > 0)) | (isempty(MA) & (m1 > 0))
%
%     General ARMA conditional mean with no regression component.
%
      [AR , MA , C , unconditionalVariance]  =  arma02(y , r1 , r2, m1 , m2);

      Flag1  =  logical(0);   % Indicate an INCOMPLETE mean specification.

   end

else

   if m1 == 0        % Check for MA terms.

      if isempty(C) | (isempty(AR) & (r1 > 0)) | isempty(regress)
%
%        General ARX conditional mean model with no MA terms. Initial
%        estimates can be generated by a simple OLS regression.
%
         yLag               =  lagmatrix(y , [1:r1]);
         yLag(isnan(yLag))  =  mean(y);

         [QQ , RR]  =  qr([ones(size(X,1),1)  yLag  X] , 0);
         b          =  RR \ (QQ'*y);
         residuals  =  y - [ones(size(X,1),1)  yLag  X]*b;

         C       =  b(1);
         AR      =  b(2:r1+1);
         MA      =  [];
         regress =  b(r1+2:end);

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

      if isempty(C) | (isempty(AR) & (r1 > 0)) | isempty(MA) | isempty(regress)
%
%        General ARMAX conditional mean model. 
%
         C       =  0; 
         AR      =  zeros(r1,1);
         MA      =  zeros(m1,1); 
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

 [K , GARCH , ARCH] = garch02(p1, p2 , q1, q2 , unconditionalVariance);

 Flag2  =  logical(0);        % Indicate an INCOMPLETE variance specification.

end
%ARCH(1)=1;
%ARCH(2)=1;
%ARCH(3)=1;
%ARCH1=ARCH(1:2);
%ARCH2=ARCH(4:end);
%ARCH=[ARCH1(:);ARCH2(:)];
%
%
%
%
%
%
%
%
%
%
%
%
%%
%

%b=0.99999800000000;
%A  =  [zeros(1,1+r1+r2+m1+m2+size(X,2)+1)  ones(1,(p1+1)*(p2+1)-1)  ones(1,(q1+1)*(q2+1)-1) zeros(1,2)];

Aeq  =  [];
beq  =  [];
infinity     =  inf;
hundred      =  100;
C;
upperBounds  =  [inf*ones(1+r1+r2+m1+m2+size(X,2),1);100000; inf*ones([(p1+1)*(p2+1)-1]+(q1+1)*(q2+1)-2,1); inf*zeros(2,1) ];
%upperBounds = [];
x0  =  [C ; AR(:) ; MA(:) ; regress(:) ; K ; GARCH(:) ; ARCH(:);[0,0]'];
[ss,sss]=size(x0);
size(x0);
%[asd,fgh]=size(x0);
%x0=zeros(asd,fgh)+.000001;
%specc=garchset;
%lowerBounds  =  [-hundred(ones(1+r1+r2+m1+m2,1)) ; -infinity(ones(size(X,2),1)) ; 1e-10 ; zeros([((p1+1)*(p2+1)-1)],1);[1-.00001;1-.00001;1-.00001];zeros([((q1+1)*(q2+1)-1)]-3,1);-1*inf*zeros(2,1)];
%if (q1==1)
 %   lowerBounds  =  [-1000*hundred(ones(1+r1+r2+m1+m2,1)) ; -infinity(ones(size(X,2),1)) ; -1000000 ; -1000*hundred([((p1+1)*(p2+1)-1)+(q1+1)*(q2+1)-2],1);-1*inf*ones(2,1)];
 %end
 lowerBounds = 1000*hundred(ones(ss,sss));
 
    
size(lowerBounds);
r1;
r2;
m1;
m2;
size(x0);
size(lowerBounds);
[coefficients] = fminsearch('Egarchlfn22'  , x0 , []  , y , r2, r1 , m1 , m2, p1 , p2, q1 , q2, X);
%[coefficients , LLF, exitFlag , output , lambda] =  fmincon('garchlfn22'  , x0 , []  , [], Aeq , beq ,lowerBounds ,upperBounds ,[], [] ,y , r1 , r2 , m1 , m2, p1 , p2, q1 , q2, []);
r2;
[LLF , G , H , e , hh] = Egarchlfn22(coefficients , y , r1 , r2 , m1 , m2 , p1 ,p2 ,q1 , q2 , X);


