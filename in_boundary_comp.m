% This subroutine computes the boundary of the set grid_points_fatten after
% fattening it.

function boundary_index=in_boundary_comp(N, grid_points_fatten)

w=1;

boundary_index=zeros(N*N,2);

for j=1:N
     for k=1:N
             if(grid_points_fatten(1,3,j,k)==0)
                 if(j>1 && k>1 && j<N && k<N)
                     if(grid_points_fatten(1,3,j-1,k-1)==1 || grid_points_fatten(1,3,j-1,k)==1 || grid_points_fatten(1,3,j-1,k+1)==1 || grid_points_fatten(1,3,j,k-1)==1 || grid_points_fatten(1,3,j,k+1)==1 || grid_points_fatten(1,3,j+1,k-1)==1 || grid_points_fatten(1,3,j+1,k)==1 || grid_points_fatten(1,3,j+1,k+1)==1)
                        boundary_index(w,:)=[j,k];
                        w=w+1;
                     end;
                
                 end;
              
             elseif(j==1 || k==1 || j==N || k==N)
                        boundary_index(w,:)=[j,k];
                        w=w+1;
             end;
     end;
end;

boundary_index=boundary_index(1:w-1,:);

