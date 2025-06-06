# Image-processing-Using-GARCH-Models-Book-Matlab-Codes
#this document, the main files related to 2D-GARCH model are introduced (which are used in chapters 1-6)

1.	Two dimensional GARCH models 

1.1.2D-GARCH



For modeling using 2D-GARCH model use garchfitt2

Used in sections 2-1, 3-1-2, 4-1,4-3, 5-2,6-2) of book 

 function [coefficients,hh,e] = garchfitt2(spec,y,X);


inputs:

o	spec: [r1,r2,m1,m2,p1,p2,q1,q2]:  r1,r2,m1,m2: degree of 2D-ARMA model (if doesn’t use ARMA model set 0 0 0 0,    p1,p2,q1,q2: degree of 2D-GARCH model usually [1 1 1 1], so  spec is usually [0 0 0 0 1 1 1 1];
o	y : 2-D data that should be modeled using 2D-GARCH
o	X=[]
•	Important Outputs:
o	Coefficients: coefficients of fitted 2D-GARCH model
o	hh: conditional standard deviation 
