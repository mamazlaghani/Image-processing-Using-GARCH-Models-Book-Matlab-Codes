function x=normalmin(m,data);
x=-prod(exp(-((data-m).^2)/2));
