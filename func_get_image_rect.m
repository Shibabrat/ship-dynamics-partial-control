function grid_points_image = func_get_image_rect(e,R,grid_points_base)

%FUNC_GET_IMAGE_RECT computes the image of the points in a rectangular
%domain stored in the 4-D array under the map f.
%Inputs:
%       e,R:    Energy of the tube, ratio of pitch/roll natural
%       frequencies
%       grid_points_base:   Initial set of points in the domain, D
%
%   Shibabrat Naik, Last modified: 02 May 2015
%   Send error and bug reports to shiba@vt.edu

    Ny = size(grid_points_base,3);
    Nvy = size(grid_points_base,4);
   
    returnMapTime = zeros(Ny*Nvy,1);   
    trajs = struct([]);
    grid_points_image = zeros(1,5,Ny,Nvy);
    
    t = 0;
    dT = 10000;       %Choose to obtain all the intersections
    
    options = odeset('RelTol',3.1e-10,'AbsTol',1e-10, ...
        'Events',@func_event_U1);
%     options = odeset('RelTol',3.1e-14,'AbsTol',1e-14, ...
%         'Events',@func_event_U1);
    
    tic
    for jj = 1:Nvy
        for ii = 1:Ny
            y = grid_points_base(1,1,ii,jj);
            vy = grid_points_base(1,2,ii,jj);
            if grid_points_base(1,3,ii,jj) == 1,
                vx = sqrt(2*(e - 0.5*(2/R^2)*(vy^2) - (y^2)));
                ic = [0,y,vx,vy];
                [tSol, ySol, te, ye, ie] = ode45(@ship_pitch_roll, ...
                    [t t+dT], ic, options);
                grid_points_image(1,1:2,ii,jj) = [ye(end,2), ye(end,4)];
                grid_points_image(1,3,ii,jj) = 1;
                returnMapTime(ii*jj,1) = te(end);
                trajs(ii*jj).time = tSol;
                trajs(:,ii*jj).y = ySol;
                disp(ii*jj); disp(te(end));
            else
                grid_points_image(1,1:2,ii,jj) = [y, vy];
                grid_points_image(1,3,ii,jj) = 1;
                returnMapTime(ii*jj,1) = 0;
                trajs(ii*jj).time = 0;
                trajs(:,ii*jj).y = 0;
                disp(ii*jj); 
            end
        end
    end
            
            
    %         plot(yout(:,1),yout(:,2),'-b', ...
    %             yout(end,1),yout(end,2),'xr');hold on
            
    %     plot(yout(:,1),yout(:,2),'-b')

    genImageTime = toc    
    save('initial_set_image','grid_points_image');    
    save(['pitch_roll_Q_',num2str(Ny),...
        '_',num2str(Nvy),'_image_v1.mat'],'-v7.3')
    
    
end
function [stateVarsDot] = ship_pitch_roll(t, stateVars)

    R = 1.6;

    stateVarsDot = zeros(4,1);
    x = stateVars(1); 
    y = stateVars(2);
    vx = stateVars(3);
    vy = stateVars(4);
    
    stateVarsDot(1) = vx;
    stateVarsDot(2) = vy;
    stateVarsDot(3) = -x + 2*(x.*y);
    stateVarsDot(4) = -R^2*y + (R^2/2)*(x.^2);
    
end
function [value, isterminal, direction] = func_event_U1(t,y)

%event_U1x defines the event of trajectory passing the plane 
%U_1 = \{(y,v_y)|x=0; v_x > 0\}

%     if abs(t) > 1e-2
%         isterminal = 1;  %terminate after waiting for a short time
%     else
%         isterminal = 0;
%     end
            
    value = y(1);
    isterminal = 1;
    direction = 1;
    
end

