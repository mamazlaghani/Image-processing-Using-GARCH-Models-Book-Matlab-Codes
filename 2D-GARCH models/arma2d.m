function  [meanAR,meanMA,meanc,meanv,AR,MA,constant,variance]=arma2d(X,m,n,R,M)
y=blkproc(X,[m n],'function1');
[a,b]=size(y);
for i=1:b
    [AR(:,i), MA(:,i), constant(:,i), variance(:,i)]=arma0(y(:,i),R,M);
end
num=0;
num1=0;
for i=1:b
    if AR(:,i)==zeros(R,1)
        num=num+1;
    end
    if MA(:,i)==zeros(M,1)
        num1=num1+1;
    end
end
meanAR=mean(AR')*(b)/(b-num)';
meanMA=mean(MA')*b/(b-num1)';
meanc=mean(constant);
meanv=mean(variance);

