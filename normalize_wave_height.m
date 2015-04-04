function normalH = normalize_wave_height(waveH)

% Normalize the wave height to obtain h as in Soliman-Thompson-1991
% Usage:    normalH = normalize_wave_height(waveH)

    ssH = 0.22;
    alpha0 = 0.73;
    waveLambda = 221.94;
    
    waveParamH = (waveH*alpha0*pi)/waveLambda;
    normalH = waveParamH/ssH;
    
end