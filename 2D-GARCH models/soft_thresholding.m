function xx=soft_thresholding(y,wv);
[C,S] = wavedec2(y,2,wv);

%obtaining subbands

cA2 = appcoef2(C,S,wv,2);
cH2 = detcoef2('h',C,S,2);
cV2 = detcoef2('v',C,S,2);
cD2 = detcoef2('d',C,S,2);
cH1 = detcoef2('h',C,S,1);
cV1 = detcoef2('v',C,S,1);
cD1 = detcoef2('d',C,S,1);
[m,n]=size(cH2);
[mm,nn]=size(cH1);

T_cH2=1.5*std2(cH2);
T_cV2=1.5*std2(cV2);
T_cD2=1.5*std2(cD2);
T_cH1=1.5*std2(cH1);
T_cV1=1.5*std2(cV1);
T_cD1=1.5*std2(cD1);

for i=1:m
    for j=1:n
        cH2_est(i,j)=sign(cH2(i,j))*max(0,abs(cH2(i,j))-T_cH2);
    end
end
for i=1:m
    for j=1:n
        cV2_est(i,j)=sign(cV2(i,j))*max(0,abs(cV2(i,j))-T_cV2);
    end
end
for i=1:m
    for j=1:n
        cD2_est(i,j)=sign(cD2(i,j))*max(0,abs(cD2(i,j))-T_cD2);
    end
end
for i=1:mm
    for j=1:nn
        cH1_est(i,j)=sign(cH1(i,j))*max(0,abs(cH1(i,j))-T_cH1);
    end
end
for i=1:mm
    for j=1:nn
        cV1_est(i,j)=sign(cV1(i,j))*max(0,abs(cV1(i,j))-T_cV1);
    end
end
for i=1:mm
    for j=1:nn
        cD1_est(i,j)=sign(cD1(i,j))*max(0,abs(cD1(i,j))-T_cD1);
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





