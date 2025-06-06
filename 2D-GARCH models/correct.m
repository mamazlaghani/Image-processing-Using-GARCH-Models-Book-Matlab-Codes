function G=RQSC_CN(h,deltat,p,q,f);
M=1/deltat;
N=1/h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sizeQ=N;

diagvec0=[5,6*ones(1,sizeQ-2),5];
Q0=(1/8)*gallery('tridiag',ones(1,size(diagvec0,2)-1),diagvec0,ones(1,size(diagvec0,2)-1));
Q0=Q0*eye(size(Q0));
size(Q0);

diagvec1=[1,zeros(1,sizeQ-2),-1];
Q1=gallery('tridiag',-1*ones(1,size(diagvec1,2)-1),diagvec1,ones(1,size(diagvec1,2)-1));
Q1=Q1*eye(size(Q1));
size(Q1);

diagvec2=[-3,-2*ones(1,sizeQ-2),-3];
Q2=gallery('tridiag',ones(1,size(diagvec2,2)-1),diagvec2,ones(1,size(diagvec2,2)-1));
Q2=Q2*eye(size(Q2));
size(Q2);

Qx=QP*Q1;
Qxx=QP*Q2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dp = diag(p);
Dq = diag(q);
Df = diag(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=[h:h:N*h];
t=[0:deltat:M*deltat];
tt=[.05:deltat:M*deltat-.05];
gama=-x.^3;           %%% u(x,0)=gama(x)

c(:,1)=pinv(Q0)*gama';
G=g(x,t);
GG=g(x,tt);
size(G);
size(GG);
size(Dp);
size(Q2);
for j=1:M
    tmpInv=inv((Q0-deltat/2*(1/h^2*Dp*Q2+1/(2*h)*Dq*Q1+Df*Q0)));%%%%%%%% in this example Dp , Dq and Df are constant but generally they can be change with  j
    c(:,j+1)=tmpInv*(Q0+deltat/2*(1/h^2*Dp*(Q2+1/12*Qxx)+1/(2*h)*Dq*(Q1-1/12*Qx)+Df*Q0))*c(:,j);
    c(:,j+1)=c(:,j+1)+tmpInv*(0.5.*(G(:,j)+G(:,j+1)));
    %c(:,j+1)=c(:,j+1)+tmpInv*GG(:,j)
end

size(c);