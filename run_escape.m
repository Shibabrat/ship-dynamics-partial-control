clear;clc;close all

omegaN = 0.62;
omegaE = 0.527;
phiRes = 101;
pPhiRes = 101;

%setting integation tolerance
OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14,'Events','on'); % high accuracy
xCapsize = [];
xSafe = [];
timeStates = ([]);

%Dimensional roll model for ship
% [phiMesh, pPhiMesh] = meshgrid(linspace(-0.88,0.88,phiRes), ...
%     linspace(-0.5456,0.5456,pPhiRes));
% tic;
% for i = 1:phiRes
%     for j = 1:pPhiRes
%         [t,x,te,xe,ie] = ode113('boat_roll_capsize',[0 1*(2*pi)/omegaE], ...
%             [phiMesh(i,j) pPhiMesh(i,j)],OPTIONS);
%                 
%         timeStates(i,j).x = x;
%         timeStates(i,j).t = t;
%         [i,j]
%         if ~isempty(xe)
% %             xe, te
%             plot(phiMesh(i,j), pPhiMesh(i,j),'xr');hold on
%             xCapsize = [xCapsize; [phiMesh(i,j) pPhiMesh(i,j)]];
%         else
%             plot(phiMesh(i,j), pPhiMesh(i,j),'xg');hold on
%             xSafe = [xSafe; [phiMesh(i,j) pPhiMesh(i,j)]];
%         end
%     end
% end
% runTime = toc;
% save -v7.3

%Non-dimensional roll model
omegaBar = omegaE/omegaN;
[phiMesh, pPhiMesh] = meshgrid(linspace(-0.88,0.88,phiRes), ...
    linspace(-0.88,0.88,pPhiRes));
tic;
for i = 1:phiRes
    for j = 1:pPhiRes
        [t,x,te,xe,ie] = ode113('boat_roll_nd',[0 1*(2*pi)/omegaBar], ...
            [phiMesh(i,j) pPhiMesh(i,j)],OPTIONS);
                
        timeStates(i,j).x = x;
        timeStates(i,j).t = t;
        [i,j]
        if ~isempty(xe)
%             xe, te
            xCapsize = [xCapsize; [phiMesh(i,j) pPhiMesh(i,j)]];
        else
            xSafe = [xSafe; [phiMesh(i,j) pPhiMesh(i,j)]];
        end
    end
end
runTime = toc;
save -v7.3

%Setting up the flag for compiling; if any changes are made to the source or header files, set this to true
% editFlag = true;        
% 
% if editFlag,
%     mex driver_integrate.c ...
%         -I/usr/local/include/ ...
%         -L/usr/local/lib/ -lgsl -lgslcblas -lm CFLAGS="\$CFLAGS -std=c99"
% end
% 
% tic;
% [tArray,xArray] = driver_integrate(0, 17*((2*pi)/omegaE), phiMesh, pPhiMesh);
% toc;
