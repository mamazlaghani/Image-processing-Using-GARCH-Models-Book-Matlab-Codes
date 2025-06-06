function y=gamma1(a,n);
x=gamrnd(a(1),a(2),1,n);
y=fminsearch('gammamin',[1 1],'',x);