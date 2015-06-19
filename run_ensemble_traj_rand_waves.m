%SCRIPT to compute the statistical relation between continuous forcing of
%rough seas and the discrete disturbance used in partial control

global H chi
phiRes = 5;
pPhiRes = 5;
[phiMesh, pPhiMesh] = meshgrid(linspace(-0.88,0.88,phiRes)', ...
    linspace(-0.52,0.52,pPhiRes)');
% phiMesh = -0.26;
% pPhiMesh = -0.44;
 
NEnsemble = 100;
omegaZ = 0.527;
numSamples = 1e2;
numWaveCycles = 1;

%=====================================================================%
timeVec = linspace(0,numWaveCycles*(2*pi/omegaZ),numSamples)';

for ii = 1:size(phiMesh,1)
    for jj = 1:size(phiMesh,2)
        H0 = 2;
        dH = 0;
        numH = 1;
        for dH=0:0.1:8
            %Paramters for sea waves 
            H = H0 + dH;
            chi = 90*(pi/180);


            trajEnsemble = struct([]);
            % trajEnsemble.t0 = [];
            % trajEnsemble.x0 = [];
            xInit = [phiMesh(ii,jj) pPhiMesh(ii,jj)];
            tic;
            for kk = 1:NEnsemble
                [tO,xO] = func_get_traj_random_waves(xInit);
                trajEnsemble(kk).t0 = tO;
                trajEnsemble(kk).x0 = xO;
            end
            ensembleRunTime = toc;

            %%Generate the point under the periodic forcing
            % OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14,'Events','on'); % high accuracy
            OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14); % high accuracy w/o event catching

            [tP,xP] = ode113('boat_roll_capsize',[0 1*(2*pi)/omegaZ],xInit,OPTIONS);

            HVec(numH) = H;
    
            save(['matlab_',num2str(H),'.mat'])
            
            chi0Vec = zeros(NEnsemble,1);
            for mm = 1:NEnsemble
                chi0Vec(mm,1) = norm(trajEnsemble(mm).x0(end,:) - xP(end,:));
            end
            avgChi0Vec(numH) = mean(chi0Vec);
            stdChi0Vec(numH) = std(chi0Vec);
            
            numH = numH + 1;
            
        end
        
        save(['ensemble_traj_',num2str(H),'_',num2str(ii),'_',...
            num2str(jj),'.mat'])
    end
end



