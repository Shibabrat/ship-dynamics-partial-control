function [omegaI, waveEnergySpect, fPhi, fTheta] = ...
    func_get_random_waves(time)

%FUNC_GET_RANDON_WAVES generates the random waves for rough seas
%
% Shibabrat Naik
%
    
    global H chi

    %Setting the parameters
    R = 1.6;
    alpha0 = 0.73;
    headingAngle = chi;
    dOmega = 0.01;
    omegaFinal = 2;
    omegaNPhi = 0.62;    %Natural roll frequency
    omegaNTheta = R*omegaNPhi;
    omegaZ = 0.527;     %Characteristic frequency
    g = 9.81;           % m/s^2
    U = 4*0.514444444 ; % m/s vessel speed, 1 knot = 0.514444444 m / s
    
    %=====================================================================%
    omegaInitial = dOmega;
    N = (omegaFinal - omegaInitial)/dOmega + 1;
    omegaI = linspace(omegaInitial, omegaFinal, N)';

    waveEnergySpect = zeros(N,1);
    for ii = 1:N
        tempRatio = (omegaZ/omegaI(ii))^4;
        waveEnergySpect(ii,1) = 0.11*(H^2)* ...
            (tempRatio/omegaI(ii))*exp(-0.44*tempRatio);
    end

    omegaIE = omegaI - (((omegaI.^2)*U)/g)*cos(headingAngle);
    epsilonI = 2*pi*rand(N,1);

    for kk = 1:length(time)
        fPhi(kk) = (omegaNPhi^2)*sin(headingAngle)*alpha0*(sqrt(2*dOmega)/g)*...
            sum( (omegaI.^2).*sqrt(waveEnergySpect).*sin(omegaIE*time(kk) + epsilonI) );

        fTheta(kk) = (omegaNTheta^2)*cos(headingAngle)*alpha0*(sqrt(2*dOmega)/g)*...
            sum( (omegaI.^2).*sqrt(waveEnergySpect).*sin(omegaIE*time(kk) + epsilonI) );
    end
        
end