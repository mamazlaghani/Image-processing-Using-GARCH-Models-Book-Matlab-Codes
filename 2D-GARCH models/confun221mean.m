function [c , ceq]=confun221mean(Parameters , y , pp , R , M , P , Q , X);
len=length(Parameters);
%c=[sum(Parameters(len-pp*(P+Q):len-pp))-1];
%X=[x(len-1);x(len)];
%all GARCH parameters
par=Parameters(2:len-2*pp);
%each column is related to one GARCH parameter related to one of the
%gaussian functions of GM
parr=reshape(par,P+Q+1,pp);
par1=parr((2:P+Q+1),:);
%deleteing a0 from summation 
prob=Parameters(end-2*pp+1:end-pp);
mean=Parameters(end-pp+1:end);
u=zeros(P+Q,pp);
for i=1:pp
    u(:,i)=ones(P+Q,1)*prob(i);
end
up=u.*par1;
upr=reshape(up,(P+Q)*pp,1);
%adding one extra constraint that summation of every of GARCH parameters
%(related to every gaussian density) must be less than one
%for i=1:pp
 %   b(i)=sum(par1(:,i));
 %end
%b=b';    
c=[sum(upr)-1];

%c=[ -1+Parameters(len)*Parameters(len-2)+Parameters(len)*Parameters(len-1)];
%c=[-1+ x(len-1)*x(len-7) + x(len-1)*x(len-6) + x(len-1)*x(len-5)+ x(len)*x(len-4) +x(len)*x(len-3)+x(len)*x(len-2)];
ceq=[prob'*mean];
