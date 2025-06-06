function x=gg(n);
k=0;
for i=2:n
    if rem(n,i)==0
        k=k+1;
    end
end
x=k;