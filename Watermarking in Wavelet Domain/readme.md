
In this document, the main files related to watermarking in wavelet domain using 2D-GARCH models are introduced 


# 1.Watermarking in wavelet domain
## 1.1	Insert watermark in wavelet domain:

### Used in section 5-2-1 of book
- 	Watermark insertion in D2,V2,H2 wavelet subbands 
	
- 	[X,y,wD,wH,wV,zD,zH,zV]=watermark_insert_garch(X,wdr,wv);

###	Inputs:
-	X:input image
-	Wdr: watermark to document ratio
-	wv: type of wavelet transform for example ‘db’4
###	Important Outputs:
-	y:watermarked image, 
-	wD: inserted watermark in D2
-	wH: inserted watermark in H2
-	wV: inserted watermark in V2
------------------------------------------------------------------------------------------------------------------------------------------


## 1.2	 watermark detector in a wavelet subband based on using 2D-GARCH model

Used in section 5-2-1 of book

function H=detector_garch_2D(y,w);
### Inputs:
-	Y: wavelet subband
-	W:watermark

 
 ### Important Outputs:
-	H : is 1 if watermark detect and 0 if watermark not detect 

_--------------------------------------------------------------------------------------------------------------------------------------
## 1.3	watermark detector in a wavelet subband based on using 2D-GARCH-GG model

Used in section 5-4-1 of book

watermark detector in a wavelet subband based on using 2D-GARCH-GG model

### H=detector_garch_generalized_Gaussian_2D(y,w);

### Inputs:
-	Y: wavelet subband
-	W:watermark

 
 ### Important Outputs:
-	H : is 1 if watermark detect and 0 if watermark not detect 



### Auxiliary files 
-	watermark detector in a wavelet subband based on using Gaussian distribution (for comparison)n 

### function H=detector_gaussian_2D(y,w);

### Inputs:
-	Y: wavelet subband
-	W:watermark

 
 ### Important Outputs:
 - H : is 1 if watermark detect and 0 if watermark not detect 

 ### watermark detector in a wavelet subband based on using Generalized Gaussian distribution (for comparison)

### H=detector_generalized_Gaussian_2D(y,w);


### Inputs:
-	Y: wavelet subband
-	W:watermark

 ### Important Outputs:
-	H : is 1 if watermark detect and 0 if watermark not detect 



Note: Other files are related to different test reported in the book for example computing ROC and …









