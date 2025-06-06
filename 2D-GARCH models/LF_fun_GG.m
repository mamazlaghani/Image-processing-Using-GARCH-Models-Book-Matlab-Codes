function LF = LF_fun_GG(y,van,h,s,p);
N=length(y);
s=0;
for i=1:N
    s=s+exp((quad('object_func',-1000,1000,[],[],y(i),van,h(i),s,p)));
end
LF=s;
    