function func_grid_generator()
    % This routine generates the grid of points needed to compute any Safe Set
    % or Asymptotic Safe Set for a flow. It needs to be given the following 
    % paramenters to work properly:
    %
    % xi,yi,xf,yf: the corners of the square where we want to find the safe set.
    % N: the resolution of the grid, being the total number of points N*N.
    %
    % The grid of points generated is grid_points_base(i,j,k,l), that is a four
    % dimensional array that will contain all the points of the grid and the
    % indexes associated with them. In order to make the latter computations 
    % easier we are going to associate with any point of the grid a pair of indexes.
    % Knowing the indexes associated with a point, for example q, with indexes 
    % 'k' and 'l', it is possible to recover its actual coordinates as:
    %
    % [x,y]=grid_points_base(1,1:2,k,l);
    %
    % That is, the actual coordinates are stored in the first and second position
    % of the vector grid_points_base(1,:,k,l). 
    %
    % To define grids of points using 4 dimensional arrays (as grid_points_base), 
    % it is very useful because we can easily define different sets using that
    % structure. We will do that using the value of the third possition of the 
    % vector grid_points_base(1,:,k,l) for a given point q with indexes 'k' and 'l' 
    % to indicate if that point is part of the set or not. That is, the value 
    % of grid_points_base(1,3,k,l). If it is '1' is part and if it is '0' it 
    % is not. We will label the sets with that structure in the following 
    % algorithms as grid_points_xxxxx, being xxxxx a very descriptive name for 
    % the kind of set that it is being represented. For example grid_points_base 
    % will be the set that we will start sculpting to get to find the safe sets.
    %
    % The other set that generates this routine is grid_points_image(i,j,k,l),
    % which is also a four dimensional array. It also contains indexes and
    % coordinates. Again, if we have a point q with indexes 'k' and 'l' the
    % coordinades of the image of q, that is f(q) will be stored in the first
    % and second possitions of grid_points_image;
    %
    % [fx,fy]=grid_points_image(1,1:2,k,l);

    % Author: Shibabrat Naik
    % Modified for tube dynamics of 2-DOF ship roll-pitch model
    
    format long

    %Obtaining the energy of tube
    energyFuncHandle = @(R,x,y,vx,vy)( 0.5*vx.^2 + (vy.^2)/R^2 + ...
                                    0.5*x.^2 + y.^2 - (x.^2).*y);
    
    leftExit = importdata('xeU1_stable_exit_00_left_e3.txt');
    indLeftExit = leftExit(:,3) > 0;
    R = 1.6;
    e = energyFuncHandle(R,leftExit(indLeftExit,1), ...
        leftExit(indLeftExit,2),leftExit(indLeftExit,3), ...
        leftExit(indLeftExit,4));
    e = sum(e)/length(e);
        
    % Here we specify the resolution of the grid, for example if we choose
    % N=6001, we will have 6001x6001=36012001 points.          
    Ny = 500; 
    Nvy = 500;
    
    %%Careful flagging and book-keeping of points on the SOS
    [ptsInPhaseSpace, grid_points_base] = ...
                            func_get_initial_set_rect(e,R,Ny,Nvy);
    grid_points_image = func_get_image_rect(e,R,grid_points_base);
    
    
    % Convert the points from ellipse to a rectangle

    save('initial_set','grid_points_base');

    save('initial_set_image','grid_points_image');
    
end

    
    
