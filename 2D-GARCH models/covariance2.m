function y=covariance2(t1,t2);
y=exp(-(abs(t1-t2)).^2);
