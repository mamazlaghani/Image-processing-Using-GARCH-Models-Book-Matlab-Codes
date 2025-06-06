function y=gammamin(a,x);
y=-prod(((x.^(a(1)-1))*((a(2)^a(1))/gamma(a(1)))).*exp(-a(2)*x));
