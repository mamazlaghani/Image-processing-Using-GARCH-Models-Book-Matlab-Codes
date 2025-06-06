function [H, pValue, testStatistic, criticalValue]=archtest2D(e,q1,q2,alpha);
[m,n]=size(e);
e1=e(1+q1:end,1+q2:end);
[M,N]=size(e1);
Z=zeros((q1+1)*(q2+1),M*N);
%size(Z)
f0=zeros(M*N,1);
f1=zeros(M*N,1);
for i=1:M
    for j=1:N
        %size(makez(e,i,j,q1,q2))
        Z(:,(i-1)*N+j)=makez(e,i,j,q1,q2);
        f0((i-1)*N+j)=[e1(i,j).^2]-1;
        f1((i-1)*N+j)=[e1(i,j)];
    end
end
%size(f0)
%size(Z)
Z=Z';
R2=f0'*Z*inv(Z'*Z)*Z'*f0/(f0'*f0);
TR2=M*N*R2;
testStatistic=TR2;
pValue   =  1 - chi2cdf(testStatistic , (q1+1)*(q2+1)-1);
criticalValue  =  chi2inv(1 - alpha , (q1+1)*(q2+1)-1 );
H  =  (alpha >= pValue);
        
