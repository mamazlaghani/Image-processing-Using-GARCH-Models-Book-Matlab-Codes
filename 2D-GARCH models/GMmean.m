function x=GMmean(p1,p2,p3,m1,m2,m3,s1,s2,s3,sn,y);
a1=p1*(exp((-(sn*(m1^2)+s1*(y^2))/(2*s1*sn))+([(m1*sn+y*s1)^2])/(2*s1*sn*(sn+s1)))/[(sn+s1)^.5]);

a2=p2*(exp((-(sn*(m2^2)+s2*(y^2))/(2*s2*sn))+([(m2*sn+y*s2)^2])/(2*s2*sn*(sn+s2)))/[(sn+s2)^.5]);
b1=[p1/((2*pi*(sn+s1))^.5)]*exp((-(y-m1)^2)/(2*(s1+sn)));
b2=[p2/((2*pi*(sn+s2))^.5)]*exp((-(y-m2)^2)/(2*(s2+sn)));
x=(a1+a2)/(b1+b2);