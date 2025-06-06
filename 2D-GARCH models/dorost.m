function [b11,b12,b21,b22,b31,b32,b41,b42,b51,b52,b61,b62,b71,b72,b81,b82]=dorost(b11,b12,b21,b22,b31,b32,b41,b42,b51,b52,b61,b62,b71,b72,b81,b82);
N=length(b11);
for i=1:N
    if b11(i)==1 & b12(i)==1
        b12(i)=0;
    end
    if b21(i)==1 & b22(i)==1
        b22(i)=0;
    end
    if b31(i)==1 & b32(i)==1
        b32(i)=0;
    end
    if b41(i)==1 & b42(i)==1
        b42(i)=0;
    end
    if b51(i)==1 & b52(i)==1
        b52(i)=0;
    end
    if b61(i)==1 & b62(i)==1
        b62(i)=0;
    end
    if b71(i)==1 & b72(i)==1
        b72(i)=0;
    end
    if b81(i)==1 & b82(i)==1
        b82(i)=0;
    end
end
