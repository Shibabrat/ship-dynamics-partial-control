global H chi

numH = 1;
NEnsemble = 100;
for dH=0:0.1:8
%Paramters for sea waves 
H = 2 + dH;
chi = 90*(pi/180);

%Parameters for time domain 
omegaZ = 0.527;
numSamples = 1e2;
numWaveCycles = 1;

%=====================================================================%
timeVec = linspace(0,numWaveCycles*(2*pi/omegaZ),numSamples)';

trajEnsemble = struct([]);
% trajEnsemble.t0 = [];
% trajEnsemble.x0 = [];
xInit = [0.1 0.25];
tic;
for kk = 1:NEnsemble
    [tO,xO] = func_get_traj_random_waves(xInit);
    plot(tO,xO(:,1),'-b')
    hold on
    trajEnsemble(kk).t0 = tO;
    trajEnsemble(kk).x0 = xO;
end
ensembleRunTime = toc;

%%Generate the point under the periodic forcing
% OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14,'Events','on'); % high accuracy
OPTIONS = odeset('RelTol',3e-8,'AbsTol',1e-8); % high accuracy w/o event catching

[tP,xP] = ode113('boat_roll_capsize',[0 1*(2*pi)/omegaZ],xInit,OPTIONS);

HVec(numH) = H;
numH = numH + 1;

save(['matlab_',num2str(H),'.mat'])
end

%% Compute the disturbance magnitude for each sample in the ensemble
ii = 1;
iiH = 2;
finalH = 10;
iidH = 0.1;
while iiH < finalH

    load(['matlab_',num2str(iiH),'.mat'])
    chi0Vec = zeros(NEnsemble,1);
    for mm = 1:NEnsemble
        chi0Vec(mm,1) = norm(trajEnsemble(mm).x0(end,:) - xP(end,:));
    end
    avgChi0Vec(ii) = mean(chi0Vec);
    stdChi0Vec(ii) = std(chi0Vec);
    HVec(ii) = iiH;
    
    ii = ii + 1;
    iiH = iiH + iidH;
end

