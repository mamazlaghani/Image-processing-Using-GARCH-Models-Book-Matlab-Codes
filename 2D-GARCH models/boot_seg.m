function x=boot_seg(X,b);
[n,m]=size(X);
for i=1:n-b+1
    for j=1:m-b+1
        x(:,:,((i-1)*(m-b+1)+j))=X(i:i+b-1,j:j+b-1);
    end
end
