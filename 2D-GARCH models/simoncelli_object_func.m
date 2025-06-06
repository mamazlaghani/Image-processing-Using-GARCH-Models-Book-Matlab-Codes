function f=simoncelli_object_func(x,y,van,s,p);
f=exp(-[((y-x).^2/(2*van))+(abs(x/s).^p)]);
%f=[1/([(2*pi*h*van).^.5]*[(s^2*gamma(3/p))/gamma(1/p)])]*[exp(-[((y-x).^2/(2*van))+(abs(x/(s*(h^.5))).^p)])];
