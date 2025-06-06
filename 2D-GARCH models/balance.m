function x=balance(a);
s=size(a);
sum=0;
x=0;
for i=1:s
    if a(i)==1
        sum=sum+1;
    end
end
if sum==s/2
    x=1;
else
    x=0;
end