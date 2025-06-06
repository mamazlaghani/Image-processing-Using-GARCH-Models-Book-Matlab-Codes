function [c , ceq]=confun22(Parameters , y , pp , R , M , P , Q , X, prob);
len=length(Parameters);
%c=[sum(Parameters(len-pp*(P+Q):len-pp))-1];
%X=[x(len-1);x(len)];
par=Parameters(2:len);
parr=reshape(par,P+Q+1,pp);
par1=parr((2:P+Q+1),:);
u=zeros(P+Q,pp);
for i=1:pp
    u(:,i)=ones(P+Q,1)*prob(i);
end
up=u.*par1;
upr=reshape(up,(P+Q)*pp,1);
c=[sum(upr)-1];

%c=[ -1+Parameters(len)*Parameters(len-2)+Parameters(len)*Parameters(len-1)];
%c=[-1+ x(len-1)*x(len-7) + x(len-1)*x(len-6) + x(len-1)*x(len-5)+ x(len)*x(len-4) +x(len)*x(len-3)+x(len)*x(len-2)];
ceq=[];
