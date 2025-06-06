
#For this section we used contourlet toolbox V2 (at the end of this document, the content of this toolbox described)

#In this document, the main files related to watermarking in contourlet domain using 2D-GARCH models are introduced 


#1.	Watermarking  in contourlet domain
##1.1	Insert watermark in contourlet domain:

Used in section 5-2-2 of book

function [X,yy,wD,zD,s8]=watermark_insert_garch(X,wdr,pfilt,dfilt);

 
 Inputs:
-	X:input image
-	Wdr: watermark to document ratio
-	pfilt, dfilt: filters for contourlet transform


 Important Outputs:
-	yy:watermarked image, 
-	wD: inserted watermark 

------------------------------------------------------------------------------------------------------------------------------------------


##1.2.	watermark detector in a contourlet subband based on using 2D-GARCH model

Used in section 5-2-2 of book

function H=detector_garch_2D(y,w);


Inputs:
-	Y: contourlet subband
-	W:watermark

 
 Important Outputs:
-	H : is 1 if watermark detect and 0 if watermark not detect 

_--------------------------------------------------------------------------------------------------------------------------------------

##1.3.	watermark detector in a contourlet subband based on using 2D-GARCH-GG model

used in section 5-4-1 of book 

watermark detector in a contourlet subband based on using 2D-GARCH-GG model

H=detector_garch_generalized_Gaussian_2D(y,w);

 Inputs:
-	Y: contourlet subband
-	W:watermark


 Important Outputs:
-	H : is 1 if watermark detect and 0 if watermark not detect 



-----------------------------------------------------------------------------------------------------------------------------------------
#Auxiliary files 
-	watermark detector in a wavelet subband based on using Gaussian distribution (for comparison)n 

function H=detector_gaussian_2D(y,w);

	Inputs:
-	Y: wavelet subband
-	W:watermark
•

 Important Outputs:
- H : is 1 if watermark detect and 0 if watermark not detect 
-------------------------------------------------------------------------------------------------------------------------------------

##watermark detector in a wavelet subband based on using Generalized Gaussian distribution (for comparison)

H=detector_generalized_Gaussian_2D(y,w);


###Inputs:
-	Y: wavelet subband
-	W:watermark

###Important Outputs:
-H : is 1 if watermark detect and 0 if watermark not detect 



###Note: Other files are related to different test reported in the book for example computing ROC and …
###
