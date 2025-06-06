function y=bin2(a);
N=bin2dec(a);
y=zeros(1,15);
for i=1:N
    y(15-i+1)=1;
end
