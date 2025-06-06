function [h,y]=likelihood(a,z,mean,h0);
len=length(z);
%h0=z(1);
h(1)=h0;
for i=2:len
    h(i)=a(1)+a(2)*(z(i-1)^2)+a(3)*h(i-1);
end
p=1;
for i=2:len   %2 0r 3
  p=p*gaussdist(z(i),mean,h(i));
end
y=p;
    

