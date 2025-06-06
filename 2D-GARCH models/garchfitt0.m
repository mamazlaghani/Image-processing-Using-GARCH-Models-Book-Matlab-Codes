function coefficients = garchfitt0(spec,y,X);
R=spec(1);
M=spec(2);
P=spec(3);
Q=spec(4);
%[AR , MA , C , unconditionalVariance]  =  arma0(y , R , M);
 %C       =  0; 
 %AR      =  zeros(R,1);
 %MA      =  zeros(M,1); 
 %regress =  zeros(size(X,2),1);
 %unconditionalVariance =2.231986820697785e-006;
%[K , GARCH , ARCH] = garch0(P , Q , unconditionalVariance);
%
regress  =  [];
[AR , MA , C , unconditionalVariance]  =  arma0(y , R , M);
[K , GARCH , ARCH] = garch0(P , Q , unconditionalVariance);
b=0.99999800000000;
A  =  [zeros(1,1+R+M+size(X,2)+1)  ones(1,P)  ones(1,Q)];
Aeq  =  [];
beq  =  [];
infinity     =  inf;
hundred      =  100;
lowerBounds  =  [-hundred(ones(1+R+M,1)) ; -infinity(ones(size(X,2),1)) ; 1e-10 ; zeros(P+Q,1)];
upperBounds  =  [];
x0  =  [C ; AR(:) ; MA(:) ; regress(:) ; K ; GARCH(:) ; ARCH(:)];
specc=garchset;
[coefficients , LLF, exitFlag , output , lambda] =  fmincon('garchllfn'  , x0 , A  , b, Aeq , beq ,lowerBounds ,upperBounds ,'garchnlc', specc.Optimization ,y , R , M , P , Q , X);


