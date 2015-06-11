%Test wave generators for regular waves with random phase and random waves
%sampled from a wave spectrum

lineMarker = '-g';
% set(0,'DefaultAxesFontSize',16,'DefaultTextFontSize',10)

global H chi
H = 4.94;
chi = 60*(pi/180);
alpha0 = 0.73;
omegaNPhi = 0.62;
lambda = 221.94;
omegaZ = 0.527;
numSamples = 1000;
numWaveCycles = 10;

%=====================================================================%
timeVec = linspace(0,numWaveCycles*(2*pi/omegaZ),numSamples)';
randPhase = 2*pi*rand(numSamples,1);

fPhiRandPhase = alpha0*omegaNPhi^2*pi*(H/lambda)*sin(omegaZ*timeVec + randPhase);
% figure(1)
% plot(timeVec,fPhiRandPhase,lineMarker)
% xlabel('$t$','interpreter','latex','fontsize',20)
% ylabel('$f_{\phi}(t)$','interpreter','latex','fontsize',20)
% title('Regular waves with random phase')

[fPhi, fTheta] = func_get_random_waves(timeVec');
figure(2)
plot(timeVec,fPhi,lineMarker)
xlabel('$t$','interpreter','latex')
ylabel('$f_{\phi}(t)$','interpreter','latex')
hold on

figure(3)
plot(timeVec,fTheta,lineMarker)
xlabel('$t$','interpreter','latex')
ylabel('$f_{\theta}(t)$','interpreter','latex')
hold on



