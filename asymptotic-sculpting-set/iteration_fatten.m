% This subroutine carries out the fatten operation of the forward 
% iteration. It only needs to perform the operation to pass as arguments, 
% the indexes of the boundary and the relative indexes of a ball of radius 
% disturbance plus control.



function grid_points_reachable=iteration_fatten(N,grid_points_iteration, boundary_index,indexes_reachable)

[rows_boundary,columns_boundary]=size(boundary_index);

[rows_indexes,columns_indexes]=size(indexes_reachable);



for j=1:rows_boundary
     for k=1:rows_indexes
             if(boundary_index(j,1)+indexes_reachable(k,1)>=1 && boundary_index(j,1)+indexes_reachable(k,1)<= N && boundary_index(j,2)+indexes_reachable(k,2)>=1 && boundary_index(j,2)+indexes_reachable(k,2)<= N)
                grid_points_iteration(1,3,boundary_index(j,1)+indexes_reachable(k,1),boundary_index(j,2)+indexes_reachable(k,2))=1;
             end;
     end;
end;

grid_points_reachable=grid_points_iteration;


