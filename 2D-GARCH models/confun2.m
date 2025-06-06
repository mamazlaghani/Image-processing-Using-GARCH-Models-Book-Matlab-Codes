function [c, ceq] = confun2(Parameters , y , R , M , P , Q , X);
len=length(Parameters);
c=[sum(Parameters(len-(P+Q)+1:len))-1];
ceq=[];