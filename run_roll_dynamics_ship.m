%% Integrate a grid of points 
phiRes = 10;
pPhiRes = 10;
[phiMesh,pPhiMesh] = meshgrid(linspace(-0.88,0.88,phiRes)',linspace(-0.52,0.52,pPhiRes)');
ic = [reshape(phiMesh,phiRes*pPhiRes,1),reshape(pPhiMesh,phiRes*pPhiRes,1)];

tic; 	
[xOut,xSafeFlag] = mex_integration([0,1*12],ic); 
runTime = toc

% plot(xOut(:,1),xOut(:,2),'.g','markersize',10)
% hold on
% plot(ic(:,1),ic(:,2),'.r')

%% Integrate a single trajectory over a long time with small increments
% tic; 
% ic = [0.0, 0];
% xTraj = ic;
% ti = 0;
% tTraj = ti;
% Te = 16*((2*pi)/0.527);
% dt = (Te - ti)/1000;
% while (ti + dt) <= Te,
% 	[xOut] = mex_integration([ti, ti + dt], ic); 
% 	xTraj = [xTraj; xOut];
%     tTraj = [tTraj; ti+dt];
% 	ic = xOut;
%     ti = ti + dt;
% end
% toc;
% plot(xTraj(:,1),xTraj(:,2),'-b')
% pause(1)
% plot(t,xTraj(:,1),'-r')

% %% Compute the Poincaré Map for the system
% numIterates = 300;
% phiRes = 20;
% pPhiRes = 10;

% xOutMap = zeros(phiRes*pPhiRes,2,numIterates);
% [phiMesh,pPhiMesh] = meshgrid(linspace(-1.3,1.3,phiRes)', ...
%     linspace(-1,1,pPhiRes)');
% ic = [reshape(phiMesh,phiRes*pPhiRes,1),reshape(pPhiMesh,phiRes*pPhiRes,1)];

% for i = 1:numIterates
%     tic; 
%     [xOut] = mex_integration([0,12],ic); 
%     runTime(i) = toc
%     ic = xOut;
%     xOutMap(:,:,i) = xOut;
% end
% save('strob_map_numPts_200_numIter_300_damp.mat')




