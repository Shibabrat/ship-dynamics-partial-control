% This subroutine computes the indexes of the closest point of the grid to 
% the point 'q' with coordinates 'x' and 'y' that are passed as arguments.


function [xo,yo]=iteration_extractor(N,xi,xf,yi,yf,x,y)

j=((x-xi)*(N-1)/(xf-xi))+1;
     
k=((y-yi)*(N-1)/(yf-yi))+1;

j1=floor(j);
j2=ceil(j);

k1=floor(k);
k2=ceil(k);

if (j-j1<0.5 && k-k1<0.5)
    xo=j1;
    yo=k1;
end;

if (j-j1>=0.5 && k-k1<0.5)
    xo=j2;
    yo=k1;
end;

if (j-j1<0.5 && k-k1>=0.5)
    xo=j1;
    yo=k2;
end;

if (j-j1>=0.5 && k-k1>=0.5)
    xo=j2;
    yo=k2;
end;
