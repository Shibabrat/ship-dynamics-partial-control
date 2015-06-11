function [tOut,xOut] = func_get_traj_random_waves(x)

%FUNC_GET_TRAJ_RANDOM_WAVES obtains trajectories for the ship roll model
%with random wave forcing. A simple linear interpolation is done for the
%generated wave during integration.
%
% Shibabrat Naik
%

    global H chi

%Paramters for sea waves 
%     H = 4.94;
%     chi = 90*(pi/180);

%Parameters for time domain 
    omegaZ = 0.527;
    numSamples = 1e2;
    numWaveCycles = 1;

%=====================================================================%
timeVec = linspace(0,numWaveCycles*(2*pi/omegaZ),numSamples)';
    
%Generate random waves
[~,~,fPhi,~] = func_get_random_waves(timeVec');
    
% OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14,'Events','on'); % high accuracy
OPTIONS = odeset('RelTol',3e-8,'AbsTol',1e-8); % high accuracy w/o event catching

%Dimensional roll model for ship
[tOut,xOut] = ode113(@(t,stateVars) ship_roll_random_waves(t, stateVars, timeVec, fPhi, flag),timeVec,x,OPTIONS);
                
end
function [stateVarsDot,termEvent, direcEvent] = ...
    ship_roll_random_waves(t, stateVars, timeVec, fPhi, flag)
    
%Roll angle of vanishing stability
    phiCritical = 0.88; %From Soliman-Thompson[1991]

%Parameters for Edith Terkol
    b1 = 0.0043;
    b2 = 0.0225;
    c1 = 0.384;
    c2 = 0.1296;
    c3 = 1.0368;
    c4 = -4.059;
    c5 = 2.4052;
%     I = 1174;
%     wM = 0;

%Interpolate from the generated sea waves
    ft = interp1(timeVec,fPhi,t);
    
% if (nargin < 3 || isempty(flag)) %without event 

    phi = stateVars(1,:); 
    pPhi = stateVars(2,:);
    
    stateVarsDot(1,:) = pPhi;
    stateVarsDot(2,:) = -b1*pPhi - b2.*abs(pPhi).*pPhi - c1*phi - ...
        c2.*abs(phi).*phi - c3.*phi.^3 -c4.*abs(phi).*phi.^3 - ...
        c5.*phi.^5 + ft;
    
% else
   
%     switch lower(flag)
%         case 'events'
%             if abs(t) > 1e-2
%                 isterminal = 1;  %terminate after waiting for a short time
%             else
%                 isterminal = 0;
%             end
            
%             isterminal = 1; %1 if the integration is to terminate at a zero of this event function
%             direction = 0;
            
            %check for crossing the horizontal line through eq. point(1,0.5)
%             stateVarsDot = stateVars(1) - phiCritical;
%             stateVarsDot = abs(stateVars(1)) - phiCritical;
%             termEvent = isterminal;
%             direcEvent = direction;            
%     end
%     
% end

end