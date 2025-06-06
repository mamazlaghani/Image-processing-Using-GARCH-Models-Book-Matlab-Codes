function [x,x1,a1,a2,a3,bA2x,bH2x,bV2x,bD2x,bH1x,bV1x,bD1x]=Edenoisingun_mean3(y,p1,p2,q1,q2,wv );
%wavelet transform of y (image+noise)
%supposing two level wavelet

[C,S] = wavedec2(y,3,wv);

%obtaining subbands

cA2 = appcoef2(C,S,wv,3);
cD3 = detcoef2('d',C,S,3);
cH3 = detcoef2('h',C,S,3);
cV3 = detcoef2('v',C,S,3);
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
[aH3,bbH3,c2]=Egarchfitt2([0 0 0 0 p1 p2 q1 q2],cH3,[]);
[aV3,bbV3,c2]=Egarchfitt2([0 0 0 0 p1 p2 q1 q2],cV3,[]);
[aD3,bbD3,c2]=Egarchfitt2([0 0 0 0 p1 p2 q1 q2],cD3,[]);

%calculating variances from standard deviation
bA2=bbA2.^2;
bH2=bbH2.^2;
bV2=bbV2.^2;
bD2=bbD2.^2;
bH1=bbH1.^2;
bV1=bbV1.^2;
bD1=bbD1.^2;
bH3=bbH3.^2;
bV3=bbV3.^2;
bD3=bbD3.^2;
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
bH3x=bH3-va;
bV3x=bV3-va;
bD3x=bD3-va;


% variance can not be negative
% setting negative element of vaiance matrices to realmin
bA2x(bA2x<=0) = realmin;
bH2x(bH2x<=0) = realmin;
bV2x(bV2x<=0) = realmin;
bD2x(bD2x<=0) = realmin;
bH1x(bH1x<=0) = realmin;
bV1x(bV1x<=0) = realmin;
bD1x(bD1x<=0) = realmin;
bH3x(bH3x<=0) = realmin;
bV3x(bV3x<=0) = realmin;
bD3x(bD3x<=0) = realmin;


[mmm,nnn]=size(cA2);
[mm1,nn1]=size(cD2);
[mm,nn]=size(cH1);
va1=va*ones(mm,nn);
va2=va*ones(mm1,nn1);
va3=va*ones(mmm,nnn);


% this stage is denoising through MMSE
cA2_est=[(bA2x./bA2).*cA2]+[aA2(1)*(va3./bA2)];
cH2_est=[(bH2x./bH2).*cH2]+[aH2(1)*(va2./bH2)];
cV2_est=[(bV2x./bV2).*cV2]+[aV2(1)*(va2./bV2)];
cD2_est=[(bD2x./bD2).*cD2]+[aD2(1)*(va2./bD2)];
cH1_est=[(bH1x./bH1).*cH1]+[aH1(1)*(va1./bH1)];
cV1_est=[(bV1x./bV1).*cV1]+[aV1(1)*(va1./bV1)];
cD1_est=[(bD1x./bD1).*cD1]+[aD1(1)*(va1./bD1)];
cH3_est=[(bH3x./bH3).*cH3]+[aH3(1)*(va3./bH3)];
cV3_est=[(bV3x./bV3).*cV3]+[aV3(1)*(va3./bV3)];
cD3_est=[(bD3x./bD3).*cD3]+[aD3(1)*(va3./bD3)];
%[mmm,nnn]=size(cA2_est);
%[mm,nn]=size(cH1_est);
%calculating new C and S
cA2r_est=reshape(cA2_est,1,mmm*nnn);
cH3r_est=reshape(cH3_est,1,mmm*nnn);
cV3r_est=reshape(cV3_est,1,mmm*nnn);
cD3r_est=reshape(cD3_est,1,mmm*nnn);
cH2r_est=reshape(cH2_est,1,mm1*nn1);
cV2r_est=reshape(cV2_est,1,mm1*nn1);
cD2r_est=reshape(cD2_est,1,mm1*nn1);
cH1r_est=reshape(cH1_est,1,mm*nn);
cV1r_est=reshape(cV1_est,1,mm*nn);
cD1r_est=reshape(cD1_est,1,mm*nn);
cA2r=reshape(cA2,1,mmm*nnn);

C_est=[cA2r_est cH3r_est cV3r_est cD3r_est cH2r_est cV2r_est cD2r_est cH1r_est cV1r_est cD1r_est];
C_est1=[cA2r cH3r_est cV3r_est cD3r_est cH2r_est cV2r_est cD2r_est cH1r_est cV1r_est cD1r_est];
%inverse wavelet
x = waverec2(C_est,S,wv);
x1= waverec2(C_est1,S,wv);
%
%
%
a3=[aA2';aH3';aV3';aD3'];
a2=[aH2';aV2';aD2'];
a1=[aH1';aV1';aD1'];
