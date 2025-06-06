# In this document, the main files related to denoising images in curvelet domain using 2D-GARCH models are described 
#Denoising in curvelet domain 

For curvelet transform we used curvelab V2 by Candes et al
As mentioned in the book we used wrapping method so, we used fdct_wrapping_matlab folder and all of our files is in this folder but all folders of curvelab V2 are available 

## Denoising using 2D-GARCH model

Denoising using 2D-GARCH model use “denoising_curvelet” 

Used in section 4-3 of book

 function x=denoising_curvelet(y,p1,p2,q1,q2,va);


Inputs:
-	y:noisy image
-	p1,p2,q1,q2: degree of 2D-GARCH model
-	va: variance of noise

 Important Outputs:
- x:denoised images 

*Note: for reducing multiplicative noise first use log transform and then above function

If not know the variance of noise use following function 

x=denoising_curvelet_unkownvar(y,p1,p2,q1,q2);

 
 Inputs:
-y:noisy image
-	p1,p2,q1,q2: degree of 2D-GARCH model
-	va: variance of noise

 Important Outputs:
-	x:denoised images 




## Denoising using 2D-GARCH-GG model

Denosing using 2D-GARCH generalized Gaussian (2D-GARCH-GG) model use “denoising_curvelet_unkownvar_GG” 

Used in section 4-4 of book
x=denoising_curvelet_unkownvar_GG(y,p1,p2,q1,q2);


Inputs:
-	y:noisy image
-	p1,p2,q1,q2: degree of 2D-GARCH models in mixture model usually 1,1,1,1

 
 Important Outputs:
- x:denoised images 

*Note: for reducing multiplicative noise first use log transform and then above function




