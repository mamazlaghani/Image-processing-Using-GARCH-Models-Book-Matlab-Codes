function x3=soft_thresholding1(y,wv);
[C,S] = wavedec2(y,2,wv);

%obtaining subbands

cA2 = appcoef2(C,S,wv,2);
cH2 = detcoef2('h',C,S,2);
cV2 = detcoef2('v',C,S,2);
cD2 = detcoef2('d',C,S,2);
cH1 = detcoef2('h',C,S,1);
cV1 = detcoef2('v',C,S,1);
cD1 = detcoef2('d',C,S,1);


std_cH2=std2(cH2);
std_cV2=std2(cV2);
std_cD2=std2(cD2);
std_cH1=std2(cH1);
std_cV1=std2(cV1);
std_cD1=std2(cD1);



ccc=[cH2(:); cV2(:); cD2(:); cD1(:) ; cH1(:); cV1(:)];
std_c=std(ccc);
[mm,nn]=size(cH1);
[mmm,nnn]=size(cA2);
%
%
%
s=0;
p=1;

%
%
aaa=.5;
T_cH2=aaa*std_c;
T_cV2=aaa*std_c;
T_cD2=aaa*std_c;
T_cH1=aaa*std_c;
T_cV1=aaa*std_c;
T_cD1=aaa*std_c;



T_cH21=1.5*std_cH2;
T_cV21=1.5*std_cV2;
T_cD21=1.5*std_cD2;
T_cH11=1.5*std_cH1;
T_cV11=1.5*std_cV1;
T_cD11=1.5*std_cD1;

%
%
%
%
for i=1:mmm
    for j=1:nnn
        cH2_est(i,j)=sign(cH2(i,j))*max(0,abs(cH2(i,j))-T_cH2);
        cH21_est(i,j)=sign(cH2(i,j))*max(0,abs(cH2(i,j))-T_cH21);
    end
end
for i=1:mmm
    for j=1:nnn
        cV2_est(i,j)=sign(cV2(i,j))*max(0,abs(cV2(i,j))-T_cV2);
        cV21_est(i,j)=sign(cV2(i,j))*max(0,abs(cV2(i,j))-T_cV21);
    end
end
for i=1:mmm
    for j=1:nnn
        cD2_est(i,j)=sign(cD2(i,j))*max(0,abs(cD2(i,j))-T_cD2);
        cD21_est(i,j)=sign(cD2(i,j))*max(0,abs(cD2(i,j))-T_cD21);
    end
end
for i=1:mm
    for j=1:nn
        cH1_est(i,j)=sign(cH1(i,j))*max(0,abs(cH1(i,j))-T_cH1);
        cH11_est(i,j)=sign(cH1(i,j))*max(0,abs(cH1(i,j))-T_cH11);
    end
end
for i=1:mm
    for j=1:nn
        cV1_est(i,j)=sign(cV1(i,j))*max(0,abs(cV1(i,j))-T_cV1);
        cV11_est(i,j)=sign(cV1(i,j))*max(0,abs(cV1(i,j))-T_cV11);
    end
end
for i=1:mm
    for j=1:nn
        cD1_est(i,j)=sign(cD1(i,j))*max(0,abs(cD1(i,j))-T_cD1);
        cD11_est(i,j)=sign(cD1(i,j))*max(0,abs(cD1(i,j))-T_cD11);
    end
end
%
%
%
cA2_est=cA2;
cA2r_est=reshape(cA2_est,1,mmm*nnn);
cH2r_est=reshape(cH2_est,1,mmm*nnn);
cV2r_est=reshape(cV2_est,1,mmm*nnn);
cD2r_est=reshape(cD2_est,1,mmm*nnn);
cH1r_est=reshape(cH1_est,1,mm*nn);
cV1r_est=reshape(cV1_est,1,mm*nn);
cD1r_est=reshape(cD1_est,1,mm*nn);






cA21_est=cA2;
cA2r_est1=reshape(cA21_est,1,mmm*nnn);
cH2r_est1=reshape(cH21_est,1,mmm*nnn);
cV2r_est1=reshape(cV21_est,1,mmm*nnn);
cD2r_est1=reshape(cD21_est,1,mmm*nnn);
cH1r_est1=reshape(cH11_est,1,mm*nn);
cV1r_est1=reshape(cV11_est,1,mm*nn);
cD1r_est1=reshape(cD11_est,1,mm*nn);





%
%
%
%
C_est=[cA2r_est cH2r_est cV2r_est cD2r_est cH1r_est cV1r_est cD1r_est];
C_est1=[cA2r_est1 cH2r_est1 cV2r_est1 cD2r_est1 cH1r_est1 cV1r_est1 cD1r_est1];
x3 = waverec2(C_est,S,wv);
x31 = waverec2(C_est1,S,wv);

