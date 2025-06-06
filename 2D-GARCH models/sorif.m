function x=sorif(a,b,w,p)
n=length(b);
a=a;
c=1;
x=zeros(2,n);
s=zeros(1,1000);
while (c~=0) 
    c=0;
        for i=1:n
         s=0;
        for j=1:n
            s=s+x(1,j)*a(i,j);
        end            
        x(2,i)=x(1,i)+(w/a(i,i))*(b(i)-s);
        if ((x(2,i)-x(1,i)>p) | (x(2,i)-x(1,i)<(-1*p)))
            c=c+1;
        end
        x(1,i)=x(2,i)
    end
end
i