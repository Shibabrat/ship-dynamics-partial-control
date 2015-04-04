function gen_sculpting_data(N)

%GEN_SCULPTING_DATA takes the input of solutions of initial grid of points
%for a set of ODEs saved in the structured array timeStates.
%Inputs:
%       timeStates: A structured array of size NxN with fields x and t
%       N: number of grid points assumed same for both variables
%

    load matlab.mat
    phiRes = N;
    pPhiRes = N;

    for i = 1:phiRes
        for j = 1:pPhiRes
                
            if isempty(find(abs(timeStates(i,j).x(:,1)) >= 0.88)) 
                %test for safe trajectory 
                safeFlag = 1;
            else
                safeFlag = 0;
            end
            grid_points_base(1,:,i,j) = [timeStates(i,j).x(1,1), ...
                timeStates(i,j).x(1,2), safeFlag, 0, 0]; 
            grid_points_image(1,:,i,j) = [timeStates(i,j).x(end,1), ...
                timeStates(i,j).x(end,2), 0, 0, 0]; 

        end
    end
    
    save('initial_set','grid_points_base');

    save('initial_set_image','grid_points_image');

end