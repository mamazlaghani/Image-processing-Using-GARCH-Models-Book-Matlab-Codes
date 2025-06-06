function [AR, MA, constant, variance] = arma02(y , r1 , r2 , m1 , m2);
if (r1+r2+m1+m2) == 0
   AR        =  [];
   MA        =  [];
   constant  =  mean(y(:));
   variance  =  var(y(:),1);
   return
end
[a,b]=size(y);
if b==1
 [AR, MA, constant, variance] = arma0(y , r1 , m1);   
elseif a==1
        [AR, MA, constant, variance] = arma0(y' , r2 , m2); 
    else 
y1=reshape(y,a*b,1);
[AR1, MA1, constant1, variance1] = arma0(y1 , r1 , m1);
yy=y';
y2=reshape(yy,a*b,1);
[AR2, MA2, constant2, variance2] = arma0(y2 , r2 , m2);
z=[AR2', zeros(1,r1*r2)]';
zz=reshape(z,r2,r1+1)';
zzz=reshape(zz,r2*(r1+1),1);
AR=[AR1;zzz];
w=[MA2', zeros(1,m1*m2)]';
ww=reshape(w,m2,m1+1)';
www=reshape(ww,m2*(m1+1),1);
MA=[MA1;www];
constant=(constant1+constant2)/2;
variance=(variance1+variance2)/2;
end