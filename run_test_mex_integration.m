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
pPhiRes = 10;
[phiMesh,pPhiMesh] = meshgrid(linspace(-0.88,0.88,phiRes)',linspace(-0.52,0.52,pPhiRes)');
ic = [reshape(phiMesh,phiRes*pPhiRes,1),reshape(pPhiMesh,phiRes*pPhiRes,1)];

tic; 
[xOut] = mex_integration([0,12],ic); 
runTime = toc;

plot(xOut(:,1),xOut(:,2),'.g')
hold on
plot(ic(:,1),ic(:,2),'.r')

%% Integrate a single trajectory over a long time with small increments
tic; 
ic = [0, 0];
xTraj = ic;
ti = 0;
tTraj = ti;
Te = 16*((2*pi)/0.527);
dt = (Te - ti)/1000;
while (ti + dt) <= Te,
	[xOut] = mex_integration([ti, ti + dt], ic); 
	xTraj = [xTraj; xOut];
    tTraj = [tTraj; ti+dt];
	ic = xOut;
    ti = ti + dt;
end
toc;
% plot(tTraj,xTraj(:,1),'-b')

%% Compute the Poincaré Map for the system
numIterates = 300;
phiRes = 20;
pPhiRes = 10;

xOutMap = zeros(phiRes*pPhiRes,2,numIterates);
[phiMesh,pPhiMesh] = meshgrid(linspace(-1.3,1.3,phiRes)', ...
    linspace(-1,1,pPhiRes)');
ic = [reshape(phiMesh,phiRes*pPhiRes,1),reshape(pPhiMesh,phiRes*pPhiRes,1)];

for i = 1:numIterates
    tic; 
    [xOut] = mex_integration([0,12],ic); 
    runTime(i) = toc
    ic = xOut;
    xOutMap(:,:,i) = xOut;
end
save('strob_map_numPts_200_numIter_300.mat')




