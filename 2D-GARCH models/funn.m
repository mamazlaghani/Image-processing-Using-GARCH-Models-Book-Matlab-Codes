function y=fun(x);
va=(.88)^2;
[a,b,c]=garchfitt2([0 0 0 0 1 1 1 1],x,[]);
b2=b.^2;
b2c=b2-va;
b2c(b2c<=0) = realmin;
y=(b2c./b2).*x;

