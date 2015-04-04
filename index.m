% This subroutine computes the relative indexes that will allow to access all the
% points that are within a distance u_0 or \xi_0 of a given point q of the grid of points.

% Also generates figures for the index of control and disturbance

function [indexes_disturbance,indexes_control]=index(N,clusterFlag,grid_points_base,disturbance,control)

format long;

n=1;

center_point=floor(N/2);

indexes=zeros(N^2,2);

for j=1:N
     for k=1:N
             p=grid_points_base(1,1:2,center_point,center_point);
             q=grid_points_base(1,1:2,j,k);
             distancia=norm(p-q);
             if(distancia<=disturbance)
                 indexes(n,:)=[j,k];
                 n=n+1;
             end;
     end;
end;


indexes_disturbance=[indexes(1:n-1,1)-center_point,indexes(1:n-1,2)-center_point];

h = figure;

if clusterFlag,
    set(h,'Visible','off')		%When display is not allowed in compute clusters
end

plot(indexes_disturbance(1:n-1,1),indexes_disturbance(1:n-1,2),'.');

title('These are the relative indexes for the disturbance')

xlabel('$x$ index','Interpreter','latex');

ylabel('$\dot{x}$ index','Interpreter','latex');

if clusterFlag,
    saveas(h,['figure_',num2str(N+1),'.fig'],'fig')
end

n=1;

indexes=zeros(N^2,2);


for j=1:N
     for k=1:N
             p=grid_points_base(1,1:2,center_point,center_point);
             q=grid_points_base(1,1:2,j,k);
             distancia=norm(p-q);
             if(distancia<=control)
                 indexes(n,:)=[j,k];
                 n=n+1;
             end;
     end;
end;

indexes_control=[indexes(1:n-1,1)-center_point,indexes(1:n-1,2)-center_point];

h = figure;

if clusterFlag,
    set(h,'Visible','off')		%When display is not allowed in compute clusters
end

plot(indexes_control(1:n-1,1),indexes_control(1:n-1,2),'.')

title('These are the relative indexes for the control')

xlabel('$x$ index','Interpreter','latex');

ylabel('$\dot{x}$ index','Interpreter','latex');


if clusterFlag,
    saveas(h,['figure_',num2str(N+2),'.fig'],'fig')
end



