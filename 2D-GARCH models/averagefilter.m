function x=averagefilter(y,n);
x=blkproc(y,[n n],'funn1');

