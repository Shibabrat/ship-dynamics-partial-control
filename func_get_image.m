function [imagePtsInPhaseSpace, grid_points_image] = ...
                    func_get_image(ptsInPhaseSpace, grid_points_base)

%
%   Shibabrat Naik, Last modified: 26 April 2015
%

   
    returnMapTime = zeros(size(grid_points_base,3),1);   
    trajs = struct([]);
    grid_points_image = zeros(1,3,size(grid_points_base,3),1);
    imagePtsInPhaseSpace = zeros(size(grid_points_base,3),4);
    
    dT = 10000;       %Choose to obtain all the intersections
    
    options = odeset('RelTol',3.1e-10,'AbsTol',1e-10, ...
        'Events',@func_event_U1);
%     options = odeset('RelTol',3.1e-14,'AbsTol',1e-14, ...
%         'Events',@func_event_U1);
    
    t = 0;
    tic
    for j = 1:size(grid_points_base,3)
        if grid_points_base(1,3,j,1) == 1,
            ic = ptsInPhaseSpace(j,:); 
            [tSol, ySol, te, ye, ie] = ode45(@ship_pitch_roll,[t t+dT], ...
                ic, options);
    %         plot(yout(:,1),yout(:,2),'-b', ...
    %             yout(end,1),yout(end,2),'xr');hold on
            grid_points_image(1,1:2,j,1) = [ye(end,2), ye(end,4)];
            imagePtsInPhaseSpace(j,:) = ye(end,:); 
            returnMapTime(j) = te(end);
            trajs(j).time = tSol;
            trajs(:,j).y = ySol;
            disp(j); disp(te(end));
        elseif grid_points_base(1,3,j,1) == 0,
            ic = ptsInPhaseSpace(j,:); 
            grid_points_image(1,1:2,j,1) = [ic(2), ic(4)];
            imagePtsInPhaseSpace(j,:) = ic; 
            returnMapTime(j) = NaN;
            trajs(j).time = NaN;
            trajs(:,j).y = ic;
            disp(j); disp(te(end));
        end
    end
%     plot(yout(:,1),yout(:,2),'-b')

    genImageTime = toc    
    save(['pitch_roll_Q_',num2str(size(grid_points_base,3)), ...
        '_image_v3.mat'],'-v7.3')
    
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
