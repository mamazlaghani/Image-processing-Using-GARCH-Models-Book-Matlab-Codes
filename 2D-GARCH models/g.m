function y=g(x,t);
N=length(x);
M=length(t);
for i=1:N
    for j=1:M
        y(i,j)=2*(x(i)^3)*t(j) - 6*(t(j)^2)*x(i) + 6*x(i) + (x(i)^3)*(t(j)^2) - (x(i)^3) ;
    end
end
