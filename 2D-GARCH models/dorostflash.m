function [b91,b92]=dorostflash(b1,b2,b3);
N=length(b1);
b91=zeros(N,1);
b92=zeros(N,1);
for i=1:N
    if b1(i)==1 & b2(i)==1 & b3(i)==1
        b91(i)=1;
        b92(i)=1;
    end
    if b1(i)==0 & b2(i)==1 & b3(i)==1
        b91(i)=1;
        b92(i)=0;
    end
    if b1(i)==0 & b2(i)==0 & b3(i)==1
        b91(i)=0;
        b92(i)=1;
    end
    if b1(i)==0 & b2(i)==0 &b3(i)==0
        b91(i)=0;
        b92(i)=0;
    end      
end
