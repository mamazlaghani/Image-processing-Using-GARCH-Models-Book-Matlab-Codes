function [h,xx,yy]=histog2d(X,n1,n2);
N=length(X);
x=X(:,1);
y=X(:,2);
xmax=max(x);
xmin=min(x);
ymax=max(y);
ymin=min(y);
cellwidx=(xmax-xmin)/n1;
cellwidy=(ymax-ymin)/n2;
for I=1:n1+1
     xx(I)=xmin+(I-1)*cellwidx;  % containing cell centers
end
for I=1:n2+1
     yy(I)=ymin+(I-1)*cellwidy;  % containing cell centers
end
size(yy)
h=zeros(n1,n2);
for i=1:N
    for j=1:n1
        for k=1:n2
            if ( X(i,1)>=xx(j) & X(i,1)<=xx(j+1))
                if (X(i,2)>=yy(k) & X(i,2)<=yy(k+1))
                    h(j,k)=h(j,k)+1;
                end
            end
        end
    end
end
 hmax=max(max(h));
   hh=h/hmax;
   %imshow(hh);           
    