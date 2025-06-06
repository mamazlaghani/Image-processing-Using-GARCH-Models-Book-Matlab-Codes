function [h,w]=freqsamp(r,m,a,phase)
b1=[1 zeros(1,m-1) -(r)^m];
b2=[1 zeros(1,m)];
[h1,w1]=freqz(b1,b2,1000);
h1=h1';
n=size(a);
yy=zeros(1,1000);
q=zeros(n,1000);
for i=1:n
       [q((i-1)*1000+1:i*1000),w2]=freqz(a(i)*[1 0 0],[1 -2*r*cos(phase(i)) r^2],1000);
        %yy=(q((i-1)*1000+1:i*1000))+yy;
   end
   for i=1:n
           yy=(q((i-1)*1000+1:i*1000))+yy;
       end
    h=(h1).*yy;
w=w1;