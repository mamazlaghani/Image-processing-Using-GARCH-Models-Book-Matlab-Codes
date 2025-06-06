function [ee , h , y] = garchGGsim2D_ggrnd(spec, row1h, col1h, row1y, col1y, n, m, s,p);
mm=length(row1h);
mrawh=mean(row1h);
nn=length(col1h);
mcolh=mean(col1h);
mcoly=mean(col1y);
mrawy=mean(row1y);
h=zeros(n,m);
y=zeros(n,m);
h(1:nn,1)=col1h; h(nn+1:n,1)=mcolh*ones(n-nn,1);
h(1,1:mm)=row1h; h(1,mm+1:m)=mrawh*ones(1,m-mm);
y(1:nn,1)=col1y; y(nn+1:m,1)=mcoly*ones(n-nn,1);
y(1,1:mm)=row1y; y(1,mm+1:m)=mrawy*ones(1,m-mm);
%h(1:nn,1)=col1h; h(nn+1:n,1)=.0000*ones(n-nn,1);
%h(1,1:mm)=row1h; h(1,mm+1:m)=.0000*ones(1,m-mm);
%y(1:nn,1)=col1y; y(nn+1:m,1)=.0000*ones(n-nn,1);
%y(1,1:mm)=row1y; y(1,mm+1:m)=.0000*ones(1,m-mm);



%e=normrnd(0,1,n,m);
e= ggrnd(0,s,p,n,m);




p1=spec(1);
p2=spec(2);
q1=spec(3);
q2=spec(4);

c       =  spec(5);       % Conditional mean constant.
%AR      =  garchget(spec , 'AR');      % Conditional mean AR coefficients.
%MA      =  garchget(spec , 'MA');      % Conditional mean MA coefficients.

a0 = spec(6); 
a10=spec(7);
a01=spec(8);
a11=spec(9);
B10=spec(10);
B01=spec(11);
B11=spec(12);% Conditional variance constant.



for i=2:n
    for j=2:m
        h(i,j)=a0 + a10*h(i-1,j) + a01*h(i,j-1) + a11*h(i-1,j-1) + B10*(y(i-1,j)^2) + B01*(y(i,j-1)^2) + B11*(y(i-1,j-1)^2);
        y(i,j)=e(i,j)*((h(i,j)^.5));
    end
end

ee=y;
y=ee+c;