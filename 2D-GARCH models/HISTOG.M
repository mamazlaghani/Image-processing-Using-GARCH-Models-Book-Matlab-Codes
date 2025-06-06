function[h,z]=histog(x,ncell);
%function[h,z]=histog(x,ncell);
%x is the data
%ncell is number of cells
%h is the histogram 
xmax=max(x);
xmin=min(x);
cellwid=(xmax-xmin)/ncell;
N=length(x);
  for I=1:ncell
     xx(I)=xmin+cellwid/2+(I-1)*cellwid;  % containing cell centers
  end
      h=zeros(1,ncell);
     for k=1:ncell
        cellmin=xx(k)-cellwid/2;
        cellmax=xx(k)+cellwid/2;
           for I=1:N
              if ( x(I)>cellmin & x(I)<cellmax)
                 h(k)=h(k)+1/(N*cellwid);
              end
           end
     end
z=xmin:(xmax-xmin)/(ncell-1):xmax;
subplot(211)
plot(z,h)
grid
subplot(212)
hist(x,length(x))

grid
figure(gcf)