% This subroutine computes the intersection between the reachable points of 
% the forward iteration and the original safe set.


function grid_points_asymptotic=intersection(N,grid_points_base,grid_points_reachable)

grid_points_asymptotic=grid_points_base;

for j=1:N
     for k=1:N
         if(grid_points_base(1,3,j,k)==1 && grid_points_reachable(1,3,j,k)==1)
                grid_points_asymptotic(1,3,j,k)=1;
         else
                grid_points_asymptotic(1,3,j,k)=0;
         end;
     end;
     
end;
