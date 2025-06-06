function z=mse(a,x,y);
n=length(x);
for i=1:n
    yy(i)=fun(a,x(i));
end
s=0;
for i=1:n
    s=s+(y(i)-yy(i))^2;
end
z=s;
