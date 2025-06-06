function x=resample1(xx,w,N);
%
%
% objective : resampling 
% 
% inputs:
% xx = states before resampling
% w = weigths of particles before resampling
% N = num of particles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% output
% x = states after resampling with same weigths
c(1)=w(1);
for i=2:N
    c(i)=c(i-1)+w(i);
end
i=1;
u1=unifrnd(0,1/N,1,1);
for j=1:N
    u(j)=u1+(j-1)/N;
    while (u(j)>c(i))
        i=i+1;
    end
    x(j)=xx(i);
end 
x=x';