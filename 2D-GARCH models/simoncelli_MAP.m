function x=simoncelli_MAP(y,s,p,van);
% y and s are vectors because s can be changed in time
N=length(y);
x=zeros(N,1);
M=length(s);
if M==1
    s=s*ones(N,1);
end
for i=1:N
    x(i)=fminsearch('simoncelli',y(i),[],s(i),p,y(i),van);
end