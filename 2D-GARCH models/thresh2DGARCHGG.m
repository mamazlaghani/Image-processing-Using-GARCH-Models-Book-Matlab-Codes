function xx=thresh2DGARCHGG(y,p1,p2,q1,q2,wv);
yn=medfilt2(y,[2 2]);
[C,S] = wavedec2(y,2,wv);
[Cn,Sn]= wavedec2(yn,2,wv);
%obtaining subbands

cA2 = appcoef2(C,S,wv,2);
cH2 = detcoef2('h',C,S,2);
cV2 = detcoef2('v',C,S,2);
cD2 = detcoef2('d',C,S,2);
cH1 = detcoef2('h',C,S,1);
cV1 = detcoef2('v',C,S,1);
cD1 = detcoef2('d',C,S,1);


cA2n = appcoef2(Cn,Sn,wv,2);
cH2n = detcoef2('h',Cn,Sn,2);
cV2n = detcoef2('v',Cn,Sn,2);
cD2n = detcoef2('d',Cn,Sn,2);
cH1n = detcoef2('h',Cn,Sn,1);
cV1n = detcoef2('v',Cn,Sn,1);
cD1n = detcoef2('d',Cn,Sn,1);



%modeling subbands with GARCH
[aA2,bbA2,c2]=garchfitt2_GG_nonoise([0 0 0 0 p1 p2 q1 q2],cA2n,[]);
[aH2,bbH2,c2]=garchfitt2_GG_nonoise([0 0 0 0 p1 p2 q1 q2],cH2n,[]);
[aV2,bbV2,c2]=garchfitt2_GG_nonoise([0 0 0 0 p1 p2 q1 q2],cV2n,[]);
[aD2,bbD2,c2]=garchfitt2_GG_nonoise([0 0 0 0 p1 p2 q1 q2],cD2n,[]);
[aH1,bbH1,c2]=garchfitt2_GG_nonoise([0 0 0 0 p1 p2 q1 q2],cH1n,[]);
[aV1,bbV1,c2]=garchfitt2_GG_nonoise([0 0 0 0 p1 p2 q1 q2],cV1n,[]);
[aD1,bbD1,c2]=garchfitt2_GG_nonoise([0 0 0 0 p1 p2 q1 q2],cD1n,[]);

NNN=length(aA2);
sA2=aA2(NNN-1);
pA2=aA2(NNN);
ssA2=sA2*bbA2;

NNN=length(aH2);
sH2=aH2(NNN-1);
pH2=aH2(NNN);
ssH2=sH2*bbH2;

NNN=length(aV2);
sV2=aV2(NNN-1);
pV2=aV2(NNN);
ssV2=sV2*bbV2;

NNN=length(aD2);
sD2=aD2(NNN-1);
pD2=aD2(NNN);
ssD2=sD2*bbD2;

NNN=length(aH1);
sH1=aH1(NNN-1);
pH1=aH1(NNN);
ssH1=sH1*bbH1;

NNN=length(aV1);
sV1=aV1(NNN-1);
pV1=aV1(NNN);
ssV1=sV1*bbV1;

NNN=length(aD1);
sD1=aD1(NNN-1);
pD1=aD1(NNN);
ssD1=sD1*bbD1;



varA2x=ssA2.^2*(gamma(3/pA2)/gamma(1/pA2));
varH2x=ssH2.^2*(gamma(3/pH2)/gamma(1/pH2));
varV2x=ssV2.^2*(gamma(3/pV2)/gamma(1/pV2));
varD2x=ssD2.^2*(gamma(3/pD2)/gamma(1/pD2));
varH1x=ssH1.^2*(gamma(3/pH1)/gamma(1/pH1));
varV1x=ssV1.^2*(gamma(3/pV1)/gamma(1/pV1));
varD1x=ssD1.^2*(gamma(3/pD1)/gamma(1/pD1));

[m,n]=size(cH2);
[mm,nn]=size(cH1);
%calculating variance of noise
[aval,dovom]=size(cH1);
yyh=abs(reshape(cH1,aval*dovom,1));
vaa=(median(yyh))/(.6745);
va=vaa^2;

varH2x(varH2x<=0) = realmin;
varV2x(varV2x<=0) = realmin;
varD2x(varD2x<=0) = realmin;
varH1x(varH1x<=0) = realmin;
varV1x(varV1x<=0) = realmin;
varD1x(varD1x<=0) = realmin;


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





