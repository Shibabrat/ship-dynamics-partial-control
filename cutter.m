% This subroutine cuts the parts of grid_points_base that do not satisfy
% the property of a safe set.

function grid_points_cutted=cutter(N,grid_points_base,grid_points_image,grid_points_shrink,xi,xf,yi,yf)

grid_points_cutted=grid_points_base;

for j=1:N
     for k=1:N
         if(grid_points_base(1,3,j,k)==1)
             q=grid_points_image(1,1:2,j,k);
             
             if(q(1)>=xi && q(1)<=xf && q(2)>=yi && q(2)<=yf)
                 [j1,j2,k1,k2]=extractor(N,xi,xf,yi,yf,q(1),q(2));
                 
                 if(grid_points_shrink(1,3,j1,k1)==1 || grid_points_shrink(1,3,j1,k2)==1 || grid_points_shrink(1,3,j2,k1)==1 || grid_points_shrink(1,3,j2,k2)==1)
                     grid_points_cutted(1,3,j,k)=1;
                 else
                     grid_points_cutted(1,3,j,k)=0;
                 end;
                 
             else
                 grid_points_cutted(1,3,j,k)=0;
             end;
             
         end;
         
     end;
     
end;


