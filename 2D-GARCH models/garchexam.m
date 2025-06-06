function y=garchexam(m,sigma,n);
%y=zeros(n,1);
y(1)=normrnd(m,sigma,1,1);
for i=2:n
    sigma1=y(i-1)^2;
    y(i)=normrnd(m,sigma1*.01,1,1);
end

