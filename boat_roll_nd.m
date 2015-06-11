function [stateVarsDot,termEvent, direcEvent] = boat_roll_nd(t, stateVars, flag)

%Roll angle of vanishing stability
    phiCritical = 0.88;

%Parameters for Edith Terkol
omegaN = 0.62;
omegaE = 0.527;
omegaBar = omegaE/omegaN;

alpha0 = 0.73;
H = 9.84;
lambda = 221.94;
b1 = 0.0043/omegaN;
b2 = 0.0225;
c2=0.1296/(omegaN^2);
c3=1.0368/(omegaN^2);
c4 = -4.059/(omegaN^2);
c5 = 2.4052/(omegaN^2);

%Disturbance due to regular seas
HBar = alpha0*pi*(H/lambda);
    
if (nargin < 3 || isempty(flag)) %without event 

    phi = stateVars(1,:); 
    pPhi = stateVars(2,:);  
    
    stateVarsDot(1,:) = pPhi;
    stateVarsDot(2,:) = - phi - c2.*abs(phi).*phi - c3.*phi.^3 ...
        - c4.*abs(phi).*phi.^3 - c5.*phi.^5 ...
        - b1*pPhi - b2.*abs(pPhi).*pPhi + HBar.*sin(omegaBar.*t);
    
else
   
    switch lower(flag)
        case 'events'
%             if abs(t) > 1e-2
%                 isterminal = 1;  %terminate after waiting for a short time
%             else
%                 isterminal = 0;
%             end
            
            isterminal = 1; %1 if the integration is to terminate at a zero of this event function
            direction = 0;
            
            %check for crossing the horizontal line through eq. point(1,0.5)
%             stateVarsDot = stateVars(1) - phiCritical;
            stateVarsDot = abs(stateVars(1)) - phiCritical;
            termEvent = isterminal;
            direcEvent = direction;            
    end
    
end

end