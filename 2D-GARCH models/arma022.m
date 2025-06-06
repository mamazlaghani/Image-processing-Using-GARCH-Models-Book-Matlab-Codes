function [AR, MA, constant, variance] = arma022(y , r1 , r2 , m1 , m2);
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
[meanAR,meanMA,meanc,meanv,AR,MA,constant,variance]=arma2d(y,max(r1,m1)+1,max(r2,m2)+1,(r1+1)*(r2+1)-1,(m1+1)*(m2+1)-1);
AR=meanAR;
MA=meanMA;
constant=meanc;
variance=meanv;
end