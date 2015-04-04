% This subroutine computes the boundary of the set grid_points_base before
% fattening it.



function boundary_indexes=in_boundary(N,grid_points_base)
w=1;

boundary_indexes=zeros(N*N,2);

for j=1:N
     for k=1:N
             if(grid_points_base(1,3,j,k)==1)
                 if(j>1 && k>1 && j<N && k<N)
                     if(grid_points_base(1,3,j-1,k-1)==0 || grid_points_base(1,3,j-1,k)==0 || grid_points_base(1,3,j-1,k+1)==0 || grid_points_base(1,3,j,k-1)==0 || grid_points_base(1,3,j,k+1)==0 || grid_points_base(1,3,j+1,k-1)==0 || grid_points_base(1,3,j+1,k)==0 || grid_points_base(1,3,j+1,k+1)==0)
                        boundary_indexes(w,:)=[j,k];
                        w=w+1;
                     end;
                     
                 elseif(j==1 && k>1 && k<N)
                     if(grid_points_base(1,3,j,k-1)==0 || grid_points_base(1,3,j,k+1)==0 || grid_points_base(1,3,j+1,k-1)==0 || grid_points_base(1,3,j+1,k)==0 || grid_points_base(1,3,j+1,k+1)==0)
                        boundary_indexes(w,:)=[j,k];
                        w=w+1;
                     end;
                     
                 elseif(k==1 && j>1 && j<N)
                     if(grid_points_base(1,3,j-1,k)==0 || grid_points_base(1,3,j-1,k+1)==0 || grid_points_base(1,3,j,k+1)==0 || grid_points_base(1,3,j+1,k)==0 || grid_points_base(1,3,j+1,k+1)==0)
                        boundary_indexes(w,:)=[j,k];
                        w=w+1;
                     end;
                 elseif(j==N && k>1 && k<N)
                     if(grid_points_base(1,3,j-1,k-1)==0 || grid_points_base(1,3,j-1,k)==0 || grid_points_base(1,3,j-1,k+1)==0 || grid_points_base(1,3,j,k-1)==0 || grid_points_base(1,3,j,k+1)==0)
                        boundary_indexes(w,:)=[j,k];
                        w=w+1;
                     end;
                 elseif(k==N && j>1 && j<N)
                     if(grid_points_base(1,3,j-1,k-1)==0 || grid_points_base(1,3,j-1,k)==0 || grid_points_base(1,3,j,k-1)==0 || grid_points_base(1,3,j+1,k-1)==0 || grid_points_base(1,3,j+1,k)==0)
                        boundary_indexes(w,:)=[j,k];
                        w=w+1;
                     end;
                 elseif(j==1 && k==1)
                     if(grid_points_base(1,3,j,k+1)==0 || grid_points_base(1,3,j+1,k)==0 || grid_points_base(1,3,j+1,k+1)==0)
                        boundary_indexes(w,:)=[j,k];
                        w=w+1;
                     end;
                        
                 elseif(j==1 && k==N)
                     if(grid_points_base(1,3,j,k-1)==0 || grid_points_base(1,3,j+1,k-1)==0 || grid_points_base(1,3,j+1,k)==0)
                        boundary_indexes(w,:)=[j,k];
                        w=w+1;
                     end;
                        
                 elseif(j==N && k==N)
                     if(grid_points_base(1,3,j-1,k)==0 || grid_points_base(1,3,j-1,k-1)==0 || grid_points_base(1,3,j,k-1)==0)
                        boundary_indexes(w,:)=[j,k];
                        w=w+1;
                     end;
                     
                 elseif(j==N && k==1)
                     if(grid_points_base(1,3,j-1,k)==0 || grid_points_base(1,3,j-1,k+1)==0 || grid_points_base(1,3,j,k+1)==0)
                        boundary_indexes(w,:)=[j,k];
                        w=w+1;
                     end;
                        
                
                 end;
             end;
     end;
end;


boundary_indexes=boundary_indexes(1:w-1,:);

