function [stateVarsDot,termEvent, direcEvent] = ...
    boat_roll_capsize(t, stateVars, flag)

    global H chi
%Roll angle of vanishing stability
    phiCritical = 0.88;
%     H = 9.84;
    
% Parameters for Edith Terkol
    b1 = 0.0043;
    b2 = 0.0225;
    c1 = 0.384;
    c2 = 0.1296;
    c3 = 1.0368;
    c4 = -4.059;
    c5 = 2.4052;
    I = 1174;
    wM = 0;

%Disturbance due to regular seas
    alpha0 = 0.73;
    omegaN = 0.62;
    omegaE = 0.527;
    lambda = 221.94;
    Mt = I*alpha0*omegaN^2*pi*(H/lambda)*(sin(chi))*sin(omegaE*t);
    
if (nargin < 3 || isempty(flag)) %without event 

%     stateVarsDot = zeros(2,1);
    phi = stateVars(1,:); 
    pPhi = stateVars(2,:);
    
    stateVarsDot(1,:) = pPhi;
    stateVarsDot(2,:) = -b1*pPhi - b2.*abs(pPhi).*pPhi - c1*phi - ...
        c2.*abs(phi).*phi - c3.*phi.^3 -c4.*abs(phi).*phi.^3 - c5.*phi.^5 + Mt/I;
    
else
   
    switch lower(flag)
        case 'events'
%             if abs(t) > 1e-2
%                 isterminal = 1;  %terminate after waiting for a short time
%             else
%                 isterminal = 0;
%             end
            
            isterminal = 0; %1 if the integration is to terminate at a zero of this event function
            direction = 0;
            
            %check for crossing the horizontal line through eq. point(1,0.5)
%             stateVarsDot = stateVars(1) - phiCritical;
            stateVarsDot = abs(stateVars(1)) - phiCritical;
            termEvent = isterminal;
            direcEvent = direction;            
    end
    
end

end