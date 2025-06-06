function y=covariance3(t1,t2);
y=exp(-abs(t2-t1)).*cos(t2-t1);