% This subroutine computes the four indexes of the a grid of points that 
% surround a point 'q' with coordinates 'x' and 'y'.

function [j1,j2,k1,k2]=extractor(N,xi,xf,yi,yf,x,y)


j=((x-xi)*(N-1)/(xf-xi))+1;
     
k=((y-yi)*(N-1)/(yf-yi))+1;

j1=floor(j);
j2=ceil(j);

k1=floor(k);
k2=ceil(k);

