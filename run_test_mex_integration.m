clear;close;clc

% Setting up the flag for compiling; if any changes are made to the source
% or header files, set this to true
    editFlag = true;        

    if editFlag,
%         mex func_get_image.c ...
%             -I/usr/local/include/ ...
%             -L/usr/local/lib/ -lgsl -lgslcblas -lm ...
%             CFLAGS="\$CFLAGS -std=c99"
        mex -v GCC='/usr/bin/gcc' mex_integration.c integration.c ...
            -I/usr/local/include/ ...
            -L/usr/local/lib/ -lgsl -lgslcblas -lm ...
            CFLAGS="\$CFLAGS -std=c99"
    end
    
%% Integrate a grid of points
phiRes = 2;
pPhiRes = 1000000;
[phiMesh,pPhiMesh] = meshgrid(linspace(-0.88,0.88,phiRes)',linspace(-0.52,0.52,pPhiRes)');
ic = [reshape(phiMesh,phiRes*pPhiRes,1),reshape(pPhiMesh,phiRes*pPhiRes,1)];

tic; 
[xOut] = mex_integration([0,12],ic); 
runTime = toc;

plot(xOut(:,1),xOut(:,2),'.k')
hold on
plot(ic(:,1),ic(:,2),'.r')

%% Integrate a single trajectory with small increments
% tic; 
% ic = [0, 0.1];
% xTraj = ic;
% ti = 0;
% Te = 48;
% dt = (Te - ti)/10000;
% while (ti + dt) <= Te,
% 	[xOut] = mex_integration([ti, ti + dt], ic); 
% 	ti = ti + dt;
% 	xTraj = [xTraj; xOut];
% 	ic = xOut;
% end
% toc;
% plot(xTraj(:,1),xTraj(:,2),'-b')
