function x=lognormalmin(m,data);
x=-prod(1/(m(2)*(2*pi)^.5)*exp(-((log(data)-m(1)).^2)/(2*(m(2)^2))));