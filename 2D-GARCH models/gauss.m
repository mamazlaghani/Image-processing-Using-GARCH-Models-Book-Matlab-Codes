function [x,del]=gauss(a,b,k)
v=clock;
n=length(a);
x=zeros(k,n);
s=zeros(n+1);
for j=2:k+1 
    for i=1:n
        for m=1:n
            if m == i
                s(m+1)=s(m);
            elseif m<i
            s(m+1)=s(m)+a(i,m)*x(j,m);
        else
              s(m+1)=s(m)+a(i,m)*x(j-1,m);
        end
        end
        x(j,i)=(1/a(i,i))*(b(i)-s(n+1));
    end
end
del=etime(clock,v);