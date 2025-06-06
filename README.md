# Image-processing-Using-GARCH-Models-Book-Matlab-Codes
This document, the main files related to 2D-GARCH model are introduced (which are used in chapters 1-6)

# 1. Two Dimensional GARCH Models

## 1.1. 2D-GARCH




For modeling using 2D-GARCH model use garchfitt2

Used in sections 2-1, 3-1-2, 4-1,4-3, 5-2,6-2) of book 

 function [coefficients,hh,e] = garchfitt2(spec,y,X);


inputs:

- spec: [r1,r2,m1,m2,p1,p2,q1,q2]:  r1,r2,m1,m2: degree of 2D-ARMA model (if doesn’t use ARMA model set 0 0 0 0, p1,p2,q1,q2: degree of 2D-GARCH model usually [1 1 1 1], so spec is usually [0 0 0 0 1 1 1 1])
- y: 2-D data that should be modeled using 2D-GARCH
- X=[]
- Important Outputs:
  - Coefficients: coefficients of fitted 2D-GARCH model
  - hh: conditional standard deviation

## 2D-GARCH-M (2D-GARCH-mixture)
For modeling using 2D-GARCH mixture model use garchfitt2_mixture

Used in sections 2-3,4-2, 4-5-1, 5-4-1, 6-4-2 of book 

function [coefficients,hh, hh1, hh2, e] = garchfitt2_mixture(spec,y,X);

inputs:
-spec: [r1,r2,m1,m2,p1,p2,q1,q2]:  r1,r2,m1,m2: degree of 2D-ARMA model (if doesn’t use ARMA model set 0 0 0 0,    p1,p2,q1,q2: degree of 2D-GARCH mixture  model usually [1 1 1 1], so  spec is usually [0 0 0 0 1 1 1 1]; number of mixture is set to 2:
-	y : 2-D data that should be modeled using 2D-GARCH
- 	X=[]

Important Outputs:
- Coefficients: coefficients of fitted 2D-GARCH models
-	hh1: conditional standard deviation of first 2D-GARCH model in the mixture 
-	hh2: conditional standard deviation of second 2D-GARCH model in the mixture 

