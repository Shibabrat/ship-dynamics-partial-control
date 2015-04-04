% This subroutine computes the relative indexes that will allow to access all the
% points that are within a distance u_0+\xi_0 of a given point q of the grid of points.

function indexes_bound=index(N,grid_points_base,disturbance,control);

format long;

n=1;

center_point=floor(N/2);

indexes=zeros(N^2,2);

for j=1:N
     
     for k=1:N
             p=grid_points_base(1,1:2,center_point,center_point);
             q=grid_points_base(1,1:2,j,k);
             distancia=norm(p-q);
             if(distancia<=disturbance+control)
                 indexes(n,:)=[j,k];
                 n=n+1;
             end;
     end;
end;

indexes=indexes(1:n-1,:);

indexes_bound=[indexes(:,1)-center_point,indexes(:,2)-center_point];

h = figure;

% set(h,'Visible','off')      %When display is not allowed in compute clusters

plot(indexes_bound(:,1),indexes_bound(:,2),'.')

title('These are the relative indexes that are reachable')

xlabel('$x$ index','Interpreter','latex');

ylabel('$\dot{x}$ index','Interpreter','latex');

% saveas(h,['figure_',num2str(N+3),'.fig'],'fig')







