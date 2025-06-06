function a=fmin(x,y);
a=fminsearch('mse',[2 12 4 5 3],'',x,y);