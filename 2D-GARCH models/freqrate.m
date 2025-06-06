function x=freqrate(y);
phi=angle(y);
f=(1/(2*pi))*diff(phi)/.01;
x=diff(f)/.01;