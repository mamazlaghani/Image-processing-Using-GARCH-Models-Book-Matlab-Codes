function [encode,decode]=dpcm(t);
y=t;
%y=wavread(t);
[pred,code,part]=dpcmopt(y,1,[-1:.1:1]);
encode=dpcmenco(y,code,part,pred);
decode=dpcmdeco(encode,code,pred);



