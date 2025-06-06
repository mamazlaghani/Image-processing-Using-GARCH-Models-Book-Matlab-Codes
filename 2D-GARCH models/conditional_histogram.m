function [h,B,xx,yy]=conditional_histogram(x,y,n1,n2);
hy=histog1(y,n2);
X=[x(:) y(:)];
[hxy,xx,yy]=histog2d(X,n1,n2);
h=zeros(n1,n2);
for i=1:n2
    if hy(i)~=0
      
    %h(:,i)=hxy(:,i);
    %else
        h(:,i)=hxy(:,i)/hy(i);
    end
end
 hmax=max(max(h));
 hh=h/hmax;
 B=hh;
 %B = imrotate(hh,90,'nearest');
 %figure;imshow(B); 
 xx=xx(1:end-1);
 yy=yy(1:end-1);
 %figure;mesh(xx,yy,h)
 %imshow(B)   