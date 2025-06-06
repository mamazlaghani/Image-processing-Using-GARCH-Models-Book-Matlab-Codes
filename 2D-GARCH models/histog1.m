function [h,z]=histog1(x,N);
xmax=max(x);
xmin=min(x);
cellwid=(xmax-xmin)/N
for I=1:N
     xx(I)=xmin+cellwid/2+(I-1)*cellwid;  % containing cell centers
end
 h=zeros(1,N);
     for k=1:N
        cellmin=xx(k)-cellwid/2;
        cellmax=xx(k)+cellwid/2;
           for I=1:length(x)
              if ( x(I)>cellmin & x(I)<cellmax)
                 h(k)=h(k)+1;
              end
           end
     end
     (xmax-xmin)/(N-1)
     z=xmin:(xmax-xmin)/(N-1):xmax;