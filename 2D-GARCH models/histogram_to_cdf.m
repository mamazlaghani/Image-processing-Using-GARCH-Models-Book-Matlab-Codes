function [cdf,xx]=histogram_to_cdf(x,h);
N=length(x);
for i=1:N-1
    xx(i)=[x(i)+x(i+1)]/2;
end
xx(N)=x(N);
for i=1:N
    cdf(i)=sum(h(1:i));
end