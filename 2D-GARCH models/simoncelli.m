function val=simoncelli(x,s,p,y,van);
%va=(s^2*gamma(3/p))/gamma(1/p);
%van=va/3;          %variance of noise is {1/3*signal variance}
val=[(y-x).^2/(2*van)]+[abs(x/s)]^p;