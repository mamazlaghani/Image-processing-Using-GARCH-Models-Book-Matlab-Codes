function h=higtog2(x,y,n1,n2);
xmax=max(x);
xmin=min(x);
ymax=max(y);
ymin=min(y);
cellwidx=(xmax-xmin)/n1;
cellwidy=(ymax-ymin)/n2;
r=[x y];
for I=1:n1
     xx(I)=xmin+cellwidx/2+(I-1)*cellwidx;  % containing cell centers
end
for I=1:n2
     yy(I)=ymin+cellwidy/2+(I-1)*cellwidy;  % containing cell centers
end
h=zeros(n1,n2);
for k=1:n1
    cellminx=xx(k)-cellwidx/2;
    cellmaxx=xx(k)+cellwidx/2;
    for I=1:length(x)
        if ( x(I)>=cellminx  & x(I)<=cellmaxx)
            
    for j=1:n2
        cellminy=yy(j)-cellwidy/2;
        cellmaxy=yy(j)-cellwidy/2;
        for I=1:length(x)
            if ( x(I)>=cellminx  & x(I)<=cellmaxx & y(I)>=cellminy & y(I)<=cellmaxy)
                    h(k,j)=h(k,j)+1;
                end
            end
        end
    end

   hmax=max(max(h));
   hh=h/hmax;
   imshow(hh);
