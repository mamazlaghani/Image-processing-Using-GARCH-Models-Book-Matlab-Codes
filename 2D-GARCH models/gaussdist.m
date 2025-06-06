function y=gaussdist(x,mean,var);
%
%
% objective : finding the value of a Gaussian pdf in point x
% 
%
%%%%%%%%%%%%55
% inputs:
% x = point 
% mean = mean of Gaussian pdf
% var = variance of Gaussian pdf
y=(1/((2*pi*var)^.5))*exp((-(x-mean)^2)/(2*var));


