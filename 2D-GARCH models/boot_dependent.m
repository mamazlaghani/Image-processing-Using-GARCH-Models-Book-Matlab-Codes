function aH=boot_dependent(X,b,k,n,m,N);
%function [aH]=boot_dependent(X,b,k,n,m,N);

%Note:
%This function is written for 2D-GARCH model. To change the model, the line
%indicated by "for 2D-GARCH model" should be changed.

%Inputs:
% "X" is the input two-dimensional data such as an image or  a subband of
% an image
%"b" is the size of subblocks use in bootstrap method if the size of image
% is M*M, [M^.5] is a fine value for b
% "n" ,"m" and "k" : It is assumed that for obtaining a two-dimensional
% image with size of approximatley, M*M, we should use n*m blocks of size (b*b) (n vertical, m horizontal). "k" is equal to n*m
%"N" : indicates the number of repeats in bootstrap method

%Outputs:
% "aH" is a matrix containing the parameter of model for all repeats

S=boot_seg(X,b);
[aa,bb,c]=size(S);
r=1:1:c;
for jj=1:N
    rr=bootrsp(r,1);
    for i=1:k
        SS(:,:,i)=S(:,:,rr(i));
    end   
    XX=boot_compose(SS,n,m);
    [aH(:,jj),bH,cH]=garchfitt2([0,0,0,0,1,1,1,1],XX,[]); % for 2D-GARCH model 
end

