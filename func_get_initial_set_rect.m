function [ptsInPhaseSpace, grid_points_base] = ...
                            func_get_initial_set_rect(e,R,Ny,Nvy)

%FUNC_GET_INITIAL_SET obtains the initial set for the pitch-roll model of
%ship dynamics using the energy manifold's intersection with the SOS, U1
%defined as
%   U_1 = \{(y,v_y) | x = 0 and v_x > 0\}
%
%   Shibabrat Naik, Last modified: 02 May 2015
%   Send error and bug reports to shiba@vt.edu

    tic
    xi=-sqrt(e);     
    yi=-R*sqrt(e);    

    xf=sqrt(e);
    yf=R*sqrt(e);
    
    y = linspace(xi,xf,Ny);
    vy = linspace(yi,yf,Nvy);
    [yRect,vyRect] = meshgrid(y,vy);
    y = reshape(yRect,Ny*Nvy,1);
    vy = reshape(vyRect,Ny*Nvy,1);
    vx = sqrt(2*(e - 0.5*(2/R^2)*(vy.^2) - (y.^2)));

    indPosRealVx = [];
    for ii = 1:length(vx)
        if isreal(vx(ii)) && vx(ii) > 0,
            indPosRealVx = [indPosRealVx; ii];
        end
    end
    y = y(indPosRealVx);
    vy = vy(indPosRealVx);
    vx = vx(indPosRealVx);
                                
    %Obtain the exits as polygon from tube dynamics computations
    leftExit = importdata('xeU1_stable_exit_00_left_e3.txt');
    indLeftExit = leftExit(:,3) > 0;

    rightExit = importdata('xeU1_stable_exit_00_right_e3.txt');
    indRightExit = rightExit(:,3) > 0;
    
    ptsInPhaseSpace = [zeros(length(y),1),y,vx,vy];
    
    grid_points_base = zeros(1,5,Ny,Nvy);    
    %Generate initial set on a rectangle circumscribing the ellipse with -1
    %flag for energetically inaccessible points, 0 for forbidden region and
    %1 for safe region.
    for jj = 1:Nvy
        for ii = 1:Ny
            grid_points_base(1,1,ii,jj) = yRect(ii,jj);
            grid_points_base(1,2,ii,jj) = vyRect(ii,jj);
            coord = [yRect(ii,jj) vyRect(ii,jj)];
            vx = sqrt(2*(e - 0.5*(2/R^2)*(vyRect(ii,jj)^2) - (yRect(ii,jj)^2)));
            if  isreal(vx) && vx > 0
                inLeftExit = inpolygon(coord(1),coord(2), ...
                    leftExit(indLeftExit,2), leftExit(indLeftExit,4));
                inRightExit = inpolygon(coord(1),coord(2), ...
                    rightExit(indRightExit,2), rightExit(indRightExit,4));
                if inLeftExit == 1 || inRightExit == 1,
                    grid_points_base(1,3,ii,jj) = 0;
%                     plot(yRect(ii,jj),vyRect(ii,jj),'.r');hold on
                else
                    grid_points_base(1,3,ii,jj) = 1;
%                     plot(yRect(ii,jj),vyRect(ii,jj),'.b');hold on
                end
            elseif grid_points_base(1,3,ii,jj) == 0
%                 plot(yRect(ii,jj),vyRect(ii,jj),'.k')
            end
        end
    end
    genPtsTime = toc
    
    save('initial_set','grid_points_base');    
    save(['pitch_roll_Q_',num2str(Ny),...
        '_',num2str(Nvy),'_base_v1.mat'])

end
