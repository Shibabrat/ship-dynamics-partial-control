function [ptsInPhaseSpace, grid_points_base] = ...
                                        func_get_initial_set(e,R,Ny,Nvy)

%FUNC_GET_INITIAL_SET obtains the initial set for the pitch-roll model of
%ship dynamics using the energy manifold's intersection with the SOS, U1
%defined as
%   U_1 = \{(y,v_y) | x = 0 and v_x > 0\}
%
%   Shibabrat Naik, Last modified: 26 April 2015
%

    tic
    xi=-sqrt(e);     
    yi=-R*sqrt(e);    

    xf=sqrt(e);
    yf=R*sqrt(e);
    
    y = linspace(xi,xf,Ny);
    vy = linspace(yi,yf,Nvy);
    [y,vy] = meshgrid(y,vy);
    y = reshape(y,Ny*Nvy,1);
    vy = reshape(vy,Ny*Nvy,1);
    vx = sqrt(2*(e - 0.5*(2/R^2)*(vy.^2) - (y.^2)));

    indPosRealVx = [];
    for ii = 1:length(vx)
        if isreal(vx(ii)) && vx(ii) >= 0,
            indPosRealVx = [indPosRealVx; ii];
        end
    end
    y = y(indPosRealVx);
    vy = vy(indPosRealVx);
    vx = vx(indPosRealVx);
    
    energyFuncHandle = @(R,x,y,vx,vy)( 0.5*vx.^2 + (vy.^2)/R^2 + ...
                                    0.5*x.^2 + y.^2 - (x.^2).*y);
                                
    %Obtain the exits as polygon from tube dynamics computations
    leftExit = importdata('xeU1_stable_exit_00_left.txt');
    indLeftExit = leftExit(:,3) > 0;

    rightExit = importdata('xeU1_stable_exit_00_right.txt');
    indRightExit = rightExit(:,3) > 0;
    
    %Flag the points which are initialized in the forbidden region
    inLeftExit = inpolygon(y,vy,leftExit(indLeftExit,2), ...
        leftExit(indLeftExit,4));
    inRightExit = inpolygon(y,vy,rightExit(indRightExit,2), ...
        rightExit(indRightExit,4));
    
%     grid_points_base = [y(~inLeftExit & ~inRightExit), ... 
%                         vx(~inLeftExit & ~inRightExit), ...
%                         vy(~inLeftExit & ~inRightExit)];

    ptsInPhaseSpace = [zeros(length(y),1),y,vx,vy];
    flagPtsInForb = ~inLeftExit & ~inRightExit;
    grid_points_base = zeros(1,3,length(y),1);
    grid_points_base(1,1,:,1) = y;
    grid_points_base(1,2,:,1) = vy;
    grid_points_base(1,3,:,1) = flagPtsInForb;
   
    genPtsTime = toc
    save(['pitch_roll_Q_',num2str(length(y)),'_base_v3.mat'])
    
    %Plotting the initial set Q and the forbidden region
    plot(y,vy,'.b')
    hold on
    plot(leftExit(indLeftExit,2),leftExit(indLeftExit,4),'.r')
    plot(rightExit(indRightExit,2),rightExit(indRightExit,4),'.r')
    plot(y(inLeftExit),vy(inLeftExit),'xg')
    plot(y(inRightExit),vy(inRightExit),'xg')
    
    
    
end
