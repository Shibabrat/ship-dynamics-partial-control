% This subroutine carries out the shrinking operation that it is needed to compute a safe set.




function grid_points_shrink=shrink(N,grid_points_fatten,boundary_index,noise_indexes)


[rows_frontera,columns_frontera]=size(boundary_index);

[rows_noise_indexes,columns_noise_indexes]=size(noise_indexes);

grid_points_shrink=grid_points_fatten;

for j=1:rows_frontera
     for k=1:rows_noise_indexes
             if(boundary_index(j,1)+noise_indexes(k,1)>=1 && boundary_index(j,1)+noise_indexes(k,1)<= N && boundary_index(j,2)+noise_indexes(k,2)>=1 && boundary_index(j,2)+noise_indexes(k,2)<= N)
                grid_points_shrink(1,3,boundary_index(j,1)+noise_indexes(k,1),boundary_index(j,2)+noise_indexes(k,2))=0;
             end;
     end;
end;
