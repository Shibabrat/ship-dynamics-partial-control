% This subroutine carries out the fatten operation that it is needed to compute a safe set.





function set_plus_control=fatten(N,grid_points_base, boundary_index,indexes_control)

[filas_frontera,columnas_frontera]=size(boundary_index);

[filas_indexes_control,columnas_indexes_control]=size(indexes_control);


for j=1:filas_frontera
     for k=1:filas_indexes_control
             if(boundary_index(j,1)+indexes_control(k,1)>=1 && boundary_index(j,1)+indexes_control(k,1)<= N && boundary_index(j,2)+indexes_control(k,2)>=1 && boundary_index(j,2)+indexes_control(k,2)<= N)
                grid_points_base(1,3,boundary_index(j,1)+indexes_control(k,1),boundary_index(j,2)+indexes_control(k,2))=1;
             end;
     end;
end;

set_plus_control=grid_points_base;