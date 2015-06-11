%SCRIPT to test admissible trajectories and simulate ship dynamics
%Admissible trajectories are based on the partial control of escape
%dynamics
global H
omegaN = 0.62;
omegaE = 0.527;
epsilonC = 2.95;
epsilonD = 3;
H = 4.94;
HBar = 0.73*pi*(H/221.94);
control = epsilonC*HBar
disturbance = epsilonD*HBar
OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14,'Events','on'); % high accuracy
numTe = 100;
Te = (2*pi)/omegaE;

%Set initial state
load initial_set.mat
ic = grid_points_base(1,1:2,100,507)

%% Now add disturbance at every period
tI = 0;
for ii = 1:numTe
    tI = (ii - 1)*Te;
    tF = ii*Te;
    %Simulate for the next encounter wave period
    [tC,xC,teC,xeC,ieC] = ode113('boat_roll_capsize',[tI tF], ...
        ic,OPTIONS);
    %Check if the state landed in the safe set 
    
    distStep = disturbance*randn(1);
    icDist = xC(end,:) + distStep*ones(1,2);
    ic = icDist;
    [tC,xC,teC,xeC,ieC] = ode113('boat_roll_capsize',[tF tF+Te], ...
        ic,OPTIONS);   
    if ieC, %capsize, so apply control step
        disp('capsize')
    else
        continue
    end   
end
                
%Fix the domain to look at the trajectories
figure(1)
plot(x(end,1),x(end,2),'or',ic(1),ic(2),'xg',x(:,1),x(:,2),'-k');hold on
xlim([-0.88 0.88])
ylim([-0.52 0.52])

figure(2)
plot(t(end),x(end,1),'or',t(1),ic(1),'xg',t,x(:,1),'-r')









