function X=boot_compose(S,n,m);
[b,b]=size(S(:,:,1));
for i=1:n
    for j=1:m
        X((i-1)*b+1:(i-1)*b+b,(j-1)*b+1:(j-1)*b+b)=S(:,:,(i-1)*m+j);
    end
end
