function LF=simoncelli_object_func_sum(para,y,van);
s=para(1);
p=para(2);
LF=1;
N=length(y);
maxy=max(y);
miny=min(y);
for i=1:N
    LF=LF*((quad('simoncelli_object_func',-500,500,[],[],y(i),van,s,p)));
end
LF=-LF;