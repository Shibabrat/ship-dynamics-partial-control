% This subroutine computes the forward iteration of the
% intermediate set. It uses the fourth element of the vector
% grid_points_asymptotic(1,4,j,k), to indicate if a particular point of the grid
% is part of the foward iteration or not. If it is '1' then is contained in
% the forward iteration, if it is '0' it is not.




function grid_points_iteration=iteration(N,grid_points_asymptotic,grid_points_image,xi,xf,yi,yf)

grid_points_iteration=grid_points_asymptotic;

grid_points_iteration(1,3,:,:)=0;

for j=1:N
     for k=1:N
         if(grid_points_asymptotic(1,3,j,k)==1)
             q=grid_points_image(1,1:2,j,k);
             
             [j1,k1]=iteration_extractor(N,xi,xf,yi,yf,q(1),q(2));
             
             grid_points_iteration(1,3,j1,k1)=1;
         end;
         
     end;
     
end;

