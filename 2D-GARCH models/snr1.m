clear all;
close all;
% it is assumed that fin=(3/N)*fs 
datab;
figure(1);
offsett=180;
N=64;  %FFT point
q=q.*boxcar(401);   % you can use the other windows. You should report the best report.
c=abs(fft(q,401));
c=c/max(c)+1e-6;

figure(2);
bar(20*log10(c(1:N/2)))
d=0;
for i=5:N/2      
 d=d+(c(i)^2);      % calculating the inband noise power
end
thd=sqrt(d)/c(4)    % c(4) represents the signal bin.
thdb=20*log10(thd) 
snrdB=-1*thdb
resolutionn_bit=snrdB/6