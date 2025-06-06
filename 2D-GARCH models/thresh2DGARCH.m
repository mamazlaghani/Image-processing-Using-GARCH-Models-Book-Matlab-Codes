function xx=thresh2DGARCH(y,p1,p2,q1,q2,wv);
[C,S] = wavedec2(y,2,wv);

%obtaining subbands

cA2 = appcoef2(C,S,wv,2);
cH2 = detcoef2('h',C,S,2);
cV2 = detcoef2('v',C,S,2);
cD2 = detcoef2('d',C,S,2);
cH1 = detcoef2('h',C,S,1);
cV1 = detcoef2('v',C,S,1);
cD1 = detcoef2('d',C,S,1);





[aA2,bbA2,c2]=garchfitt2([0 0 0 0 p1 p2 q1 q2],cA2,[]);
[aH2,bbH2,c2]=garchfitt2([0 0 0 0 p1 p2 q1 q2],cH2,[]);
[aV2,bbV2,c2]=garchfitt2([0 0 0 0 p1 p2 q1 q2],cV2,[]);
[aD2,bbD2,c2]=garchfitt2([0 0 0 0 p1 p2 q1 q2],cD2,[]);
[aH1,bbH1,c2]=garchfitt2([0 0 0 0 p1 p2 q1 q2],cH1,[]);
[aV1,bbV1,c2]=garchfitt2([0 0 0 0 p1 p2 q1 q2],cV1,[]);
[aD1,bbD1,c2]=garchfitt2([0 0 0 0 p1 p2 q1 q2],cD1,[]);


%calculating variances from standard deviation
bA2=bbA2.^2;
bH2=bbH2.^2;
bV2=bbV2.^2;
bD2=bbD2.^2;
bH1=bbH1.^2;
bV1=bbV1.^2;
bD1=bbD1.^2;


[m,n]=size(cH2);
[mm,nn]=size(cH1);
%calculating variance of noise
[aval,dovom]=size(cH1);
yyh=abs(reshape(cH1,aval*dovom,1));
vaa=(median(yyh))/(.6745);
va=vaa^2;

bA2x=bA2-va;
bH2x=bH2-va;
bV2x=bV2-va;
bD2x=bD2-va;
bH1x=bH1-va;
bV1x=bV1-va;
bD1x=bD1-va;

varH2x=bH2x;
varV2x=bV2x;
varD2x=bD2x;
varH1x=bH1x;
varV1x=bV1x;
varD1x=bD1x;

varH2x(varH2x<=0) = realmin;
varV2x(varV2x<=0) = realmin;
varD2x(varD2x<=0) = realmin;
varH1x(varH1x<=0) = realmin;
varV1x(varV1x<=0) = realmin;
varD1x(varD1x<=0) = realmin;
min(min(varH2x))

T_cH2=va./((varH2x).^.5);
T_cV2=va./((varV2x).^.5);
T_cD2=va./((varD2x).^.5);
T_cH1=va./((varH1x).^.5);
T_cV1=va./((varV1x).^.5);
T_cD1=va./((varD1x).^.5);
size(T_cH2)
size(T_cD1)

for i=1:m
    for j=1:n
        cH2_est(i,j)=sign(cH2(i,j))*max(0,abs(cH2(i,j))-T_cH2(i,j));
    end
end
for i=1:m
    for j=1:n
        cV2_est(i,j)=sign(cV2(i,j))*max(0,abs(cV2(i,j))-T_cV2(i,j));
    end
end
for i=1:m
    for j=1:n
        cD2_est(i,j)=sign(cD2(i,j))*max(0,abs(cD2(i,j))-T_cD2(i,j));
    end
end
for i=1:mm
    for j=1:nn
        cH1_est(i,j)=sign(cH1(i,j))*max(0,abs(cH1(i,j))-T_cH1(i,j));
    end
end
for i=1:mm
    for j=1:nn
        cV1_est(i,j)=sign(cV1(i,j))*max(0,abs(cV1(i,j))-T_cV1(i,j));
    end
end
for i=1:mm
    for j=1:nn
        cD1_est(i,j)=sign(cD1(i,j))*max(0,abs(cD1(i,j))-T_cD1(i,j));
    end
end

%
%
%
cA2_est=cA2;
cA2r_est=reshape(cA2_est,1,m*n);
cH2r_est=reshape(cH2_est,1,m*n);
cV2r_est=reshape(cV2_est,1,m*n);
cD2r_est=reshape(cD2_est,1,m*n);
cH1r_est=reshape(cH1_est,1,mm*nn);
cV1r_est=reshape(cV1_est,1,mm*nn);
cD1r_est=reshape(cD1_est,1,mm*nn);



C_est=[cA2r_est cH2r_est cV2r_est cD2r_est cH1r_est cV1r_est cD1r_est];
xx = waverec2(C_est,S,wv);





