function [lnx,XX]=boot_dependent1(X,b,k,n,m,N);
S=boot_seg(X,b);
[aa,bb,c]=size(S);
r=1:1:c;
for jj=1:N
    rr=bootrsp(r,1);
    for i=1:k
        SS(:,:,i)=S(:,:,rr(i));
    end   
    XX=boot_compose(SS,n,m);
    [lnx(:,:,jj),x1,a1,a2,bA2x,bH2x,bV2x,bD2x,bH1x,bV1x,bD1x]=denoisingun(XX,1,1,1,1,'db4' );
end