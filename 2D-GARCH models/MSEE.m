function mse=MSEE(a,b);
s=0;
[m,n]=size(a);
for i=1:m
    for j=1:n
        s=s+([(a(i,j)-b(i,j))^2]/(m*n));
    end
end
%mse=s/(m*n);
mse=s;
        