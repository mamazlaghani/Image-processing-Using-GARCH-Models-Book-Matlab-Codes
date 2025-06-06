function [x,x1,a1,a2,bA2x,bH2x,bV2x,bD2x,bH1x,bV1x,bD1x]=Edenoisingun(y,p1,p2,q1,q2,wv );
%wavelet transform of y (image+noise)
%supposing two level wavelet

[C,S] = wavedec2(y,2,wv);

%obtaining subbands

cA2 = appcoef2(C,S,wv,2);
cH2 = detcoef2('h',C,S,2);
cV2 = detcoef2('v',C,S,2);
cD2 = detcoef2('d',C,S,2);
cH1 = detcoef2('h',C,S,1);
cV1 = detcoef2('v',C,S,1);
cD1 = detcoef2('d',C,S,1);


%modeling subbands with GARCH
[aA2,bbA2,c2]=Egarchfitt2([0 0 0 0 p1 p2 q1 q2],cA2,[]);
[aH2,bbH2,c2]=Egarchfitt2([0 0 0 0 p1 p2 q1 q2],cH2,[]);
[aV2,bbV2,c2]=Egarchfitt2([0 0 0 0 p1 p2 q1 q2],cV2,[]);
[aD2,bbD2,c2]=Egarchfitt2([0 0 0 0 p1 p2 q1 q2],cD2,[]);
[aH1,bbH1,c2]=Egarchfitt2([0 0 0 0 p1 p2 q1 q2],cH1,[]);
[aV1,bbV1,c2]=Egarchfitt2([0 0 0 0 p1 p2 q1 q2],cV1,[]);
[aD1,bbD1,c2]=Egarchfitt2([0 0 0 0 p1 p2 q1 q2],cD1,[]);


%calculating variances from standard deviation
bA2=bbA2.^2;
bH2=bbH2.^2;
bV2=bbV2.^2;
bD2=bbD2.^2;
bH1=bbH1.^2;
bV1=bbV1.^2;
bD1=bbD1.^2;
%
%
%


%calculating variance of noise
[aval,dovom]=size(cH1);
yyh=abs(reshape(cH1,aval*dovom,1));
vaa=(median(yyh))/(.6745);
va=vaa^2;



%
%
%calculating variance of original image by subtracting variance of noise
%from variance of noisy image
bA2x=bA2-va;
bH2x=bH2-va;
bV2x=bV2-va;
bD2x=bD2-va;
bH1x=bH1-va;
bV1x=bV1-va;
bD1x=bD1-va;


% variance can not be negative
% setting negative element of vaiance matrices to realmin
bA2x(bA2x<=0) = realmin;
bH2x(bH2x<=0) = realmin;
bV2x(bV2x<=0) = realmin;
bD2x(bD2x<=0) = realmin;
bH1x(bH1x<=0) = realmin;
bV1x(bV1x<=0) = realmin;
bD1x(bD1x<=0) = realmin;


% this stage is denoising through MMSE
cA2_est=[(bA2x./bA2).*cA2];
cH2_est=[(bH2x./bH2).*cH2];
cV2_est=(bV2x./bV2).*cV2;
cD2_est=(bD2x./bD2).*cD2;
cH1_est=(bH1x./bH1).*cH1;
cV1_est=(bV1x./bV1).*cV1;
cD1_est=(bD1x./bD1).*cD1;
[mmm,nnn]=size(cA2_est);
[mm,nn]=size(cH1_est);
%calculating new C and S
cA2r_est=reshape(cA2_est,1,mmm*nnn);
cH2r_est=reshape(cH2_est,1,mmm*nnn);
cV2r_est=reshape(cV2_est,1,mmm*nnn);
cD2r_est=reshape(cD2_est,1,mmm*nnn);
cH1r_est=reshape(cH1_est,1,mm*nn);
cV1r_est=reshape(cV1_est,1,mm*nn);
cD1r_est=reshape(cD1_est,1,mm*nn);
cA2r=reshape(cA2,1,mmm*nnn);

C_est=[cA2r_est cH2r_est cV2r_est cD2r_est cH1r_est cV1r_est cD1r_est];
C_est1=[cA2r cH2r_est cV2r_est cD2r_est cH1r_est cV1r_est cD1r_est];
%inverse wavelet
x = waverec2(C_est,S,wv);
x1= waverec2(C_est1,S,wv);
%
%
%
a2=[aA2';aH2';aV2';aD2'];
a1=[aH1';aV1';aD1'];
