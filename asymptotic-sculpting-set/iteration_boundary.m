% This subroutine computes the boundary of the forward iteration of the  
% set sored in in grid_points_iteration. It stores the indexes of the 
% elements of the boundary in boundary_index.




function boundary_index=iteration_boundary(N,grid_points_iteration)
w=1;

boundary_index=zeros(N*N,2);

for j=1:N
     for k=1:N
             if(grid_points_iteration(1,3,j,k)==1)
                 if(j>1 && k>1 && j<N && k<N)
                     if(grid_points_iteration(1,3,j-1,k-1)==0 || grid_points_iteration(1,3,j-1,k)==0 || grid_points_iteration(1,3,j-1,k+1)==0 || grid_points_iteration(1,3,j,k-1)==0 || grid_points_iteration(1,3,j,k+1)==0 || grid_points_iteration(1,3,j+1,k-1)==0 || grid_points_iteration(1,3,j+1,k)==0 || grid_points_iteration(1,3,j+1,k+1)==0)
                        boundary_index(w,:)=[j,k];
                        w=w+1;
                     end;
                     
                 elseif(j==1 && k>1 && k<N)
                     if(grid_points_iteration(1,3,j,k-1)==0 || grid_points_iteration(1,3,j,k+1)==0 || grid_points_iteration(1,3,j+1,k-1)==0 || grid_points_iteration(1,3,j+1,k)==0 || grid_points_iteration(1,3,j+1,k+1)==0)
                        boundary_index(w,:)=[j,k];
                        w=w+1;
                     end;
                     
                 elseif(k==1 && j>1 && j<N)
                     if(grid_points_iteration(1,3,j-1,k)==0 || grid_points_iteration(1,3,j-1,k+1)==0 || grid_points_iteration(1,3,j,k+1)==0 || grid_points_iteration(1,3,j+1,k)==0 || grid_points_iteration(1,3,j+1,k+1)==0)
                        boundary_index(w,:)=[j,k];
                        w=w+1;
                     end;
                 elseif(j==N && k>1 && k<N)
                     if(grid_points_iteration(1,3,j-1,k-1)==0 || grid_points_iteration(1,3,j-1,k)==0 || grid_points_iteration(1,3,j-1,k+1)==0 || grid_points_iteration(1,3,j,k-1)==0 || grid_points_iteration(1,3,j,k+1)==0)
                        boundary_index(w,:)=[j,k];
                        w=w+1;
                     end;
                 elseif(k==N && j>1 && j<N)
                     if(grid_points_iteration(1,3,j-1,k-1)==0 || grid_points_iteration(1,3,j-1,k)==0 || grid_points_iteration(1,3,j,k-1)==0 || grid_points_iteration(1,3,j+1,k-1)==0 || grid_points_iteration(1,3,j+1,k)==0)
                        boundary_index(w,:)=[j,k];
                        w=w+1;
                     end;
                 elseif(j==1 && k==1)
                     if(grid_points_iteration(1,3,j,k+1)==0 || grid_points_iteration(1,3,j+1,k)==0 || grid_points_iteration(1,3,j+1,k+1)==0)
                        boundary_index(w,:)=[j,k];
                        w=w+1;
                     end;
                        
                 elseif(j==1 && k==N)
                     if(grid_points_iteration(1,3,j,k-1)==0 || grid_points_iteration(1,3,j+1,k-1)==0 || grid_points_iteration(1,3,j+1,k)==0)
                        boundary_index(w,:)=[j,k];
                        w=w+1;
                     end;
                        
                 elseif(j==N && k==N)
                     if(grid_points_iteration(1,3,j-1,k)==0 || grid_points_iteration(1,3,j-1,k-1)==0 || grid_points_iteration(1,3,j,k-1)==0)
                        boundary_index(w,:)=[j,k];
                        w=w+1;
                     end;
                     
                 elseif(j==N && k==1)
                     if(grid_points_iteration(1,3,j-1,k)==0 || grid_points_iteration(1,3,j-1,k+1)==0 || grid_points_iteration(1,3,j,k+1)==0)
                        boundary_index(w,:)=[j,k];
                        w=w+1;
                     end;
                        
                
                 end;
             end;
     end;
end;

boundary_index=boundary_index(1:w-1,:);



