function [y,count]=thresh_garch(x,h,TH);
[M,N]=size(x);
for i=1:M
    for j=1:N
        if h(i,j)<TH
            x(i,j)=0;
        end
    end
end
y=x;
count=0;
for i=1:M
    for j=1:N
        if y(i,j)==0
            count=count+1;
        end
    end
end