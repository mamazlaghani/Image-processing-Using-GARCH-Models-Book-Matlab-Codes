function [y,t]=fmm(t0,v,Ac);
t=[0:.01:50];
q=normrnd(0,v,1,1);
phi=(t-t0).^2;
fc=(1000000000)+(100)*q;
y=Ac*cos((2*pi*fc*t));
