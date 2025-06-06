function y=gaussianmixture(x,p,m,h);
N=length(m);
M=length(x);
s=0;
%
%
%
%
for j=1:M
    s=0;
    for i=1:N
        s=s+(p(i)*(1/((2*pi*h(i))^.5))*exp(-[(x(j)-m(i))^2]/(2*h(i))));
    end
    y(j)=s;
end
    
        
        
  
        
    
