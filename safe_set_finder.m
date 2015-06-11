function safe_set_finder(N,clusterFlag)
% This routine computes the safe set (if it exists) of a flow. To make it
% work properly it has to be given the following parameters:
% control: the value of the bound of the control.
% disturbance: the bound of the value of the disturbance.
% xi,yi,xf,xf: the corners of the square where we want to find the safe set.
% N: the resolution of the grid, being the total number of points N*N.
%
% It uses these three grid of points:
%
% -The array grid_points_base that will contain the initial set that is
% going to be sculpted and its successive cuts.
% -The array grid_points_fatten that will host the fattened resulting set
% during each iteration.
% -The array grid_points_shrink that will host the shrinked resulting set
% during each iteration.

% Adapted and modified by Shibabrat Naik, Shane Ross
% Date: 22 Nov 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CONFIGURATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long;

% control=0.0475;    %Here we specify the bound of the control that we are going to use.
% disturbance=0.08;    %Here we specify the bound of the disturbance that we are going to use.

% global epsilonC
% % epsilonC = 2.95;
% epsilonD = 3;
% H = 4.94;
% HBar = 0.73*pi*(H/221.94);
% control = epsilonC*HBar
% disturbance = epsilonD*HBar

% omegaN = 0.62;
% 
% xi=-0.88;    %Here are the corners of the box where we are going to 
% yi=-0.52;    %compute the Safe Sets and  the Asymptotic Safe Set.
% 
% xf=0.88;
% yf=0.52;

%Obtaining the energy of tube
energyFuncHandle = @(R,x,y,vx,vy)( 0.5*vx.^2 + (vy.^2)/R^2 + ...
                                0.5*x.^2 + y.^2 - (x.^2).*y);

leftExit = importdata('xeU1_stable_exit_00_left_e1.txt');
indLeftExit = leftExit(:,3) > 0;
R = 1.6;
e = energyFuncHandle(R,leftExit(indLeftExit,1), ...
    leftExit(indLeftExit,2),leftExit(indLeftExit,3), ...
    leftExit(indLeftExit,4));
e = sum(e)/length(e);

R = 1.6;
xi=-sqrt(e);     
yi=-R*sqrt(e);    

xf=sqrt(e);
yf=R*sqrt(e);
global control disturbance
% control=0.041;    %Here we specify the bound of the control that we are going to use.
% disturbance=0.1;    %Here we specify the bound of the disturbance that we are going to use.


safeRatio = control/disturbance
if nargin == 0 	%Default case
	N=3001;    	%Here we specify the resolution of the grid, for example if we choose
           		%N=6001, we will have 6001x6001=36012001 points.
end

% save(['safe_set_sculpting_params_',num2str(epsilonC),'.mat'])
save(['safe_set_sculpting_params_',num2str(e),'_',num2str(R),'_',num2str(control),'.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CONFIGURATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%====LCC Cluster specific===============%
if clusterFlag,
    load('~/partial-control-roll-capsize/initial_set.mat');
    load('~/partial-control-roll-capsize/initial_set_image.mat');
else
    load initial_set; %With this we load the grid of points generated with integrator and adapted to matlab format with adaptador.
    load initial_set_image; %With this we load the image of the grid of points generated with integrator and adapted to matlab format with adaptador.
end
%=======================================%

[indexes_disturbance,indexes_control]=index(N,clusterFlag,grid_points_base,disturbance,control);    %Here we compute the relative indexes needed to carry out the fatten operation and the shrink operation.

previous_num_points=-1;

fig_points=zeros(N*N,2);    %In this array we will store the points that remain after each cut in order to plot them.


for i=1:100

tic;
    
fprintf('\n\nThis is the cut number %d\n\n',i);

disp('Step 1');
boundary_index=in_boundary(N,grid_points_base);    %Here we compute the boundary of the set that remains after every cut and is stored in grid_points_base.
disp('Step 2');
grid_points_fatten=fatten(N,grid_points_base,boundary_index,indexes_control);    %Here we carry out the fatten operation.
disp('Step 3');
boundary_index_comp=in_boundary_comp(N,grid_points_fatten);    %Here we compute the boundary of the set that remains after the fatten operation.
disp('Step 4');
grid_points_shrink=shrink(N,grid_points_fatten,boundary_index_comp,indexes_disturbance);    %Here we carry out the shrink operation.
disp('Step 5');
grid_points_base=cutter(N,grid_points_base,grid_points_image,grid_points_shrink,xi,xf,yi,yf);    %Here we remove the points that do not satisfy the safe set property.



current_num_points=1;    %We use this variable to count the numbers of points that remain after every cut.


for j=1:N    %Here we count the number of points that remain after every cut.
     for k=1:N
             if(grid_points_base(1,3,j,k)==1)
                 fig_points(current_num_points,:)=grid_points_base(1,1:2,j,k);
                 current_num_points=current_num_points+1;
             end;
     end;
end;

fig_points=fig_points(1:current_num_points-1,:);


mem_points(i)=current_num_points    %In this array are stored the succesive number of points that remain after every cut.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PRINT ZONE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=figure;

if clusterFlag,
    set(h,'Visible','off')		%When display is not allowed in compute clusters
end

% plot(fig_points(:,1),fig_points(:,2),'.','MarkerSize',1)
plot(fig_points(:,1),fig_points(:,2),'.b')

iteration_number=num2str(i);

title(['Iteration number ' iteration_number ' of the Sculpting Algorithm'])

% rectangle('Position',[xi+4*disturbance,yi+8*disturbance,2*disturbance,2*disturbance],'Curvature',[1,1],'EdgeColor',[0 0 0],'FaceColor','g');
% rectangle('Position',[xi+4*disturbance+disturbance-control,yi+4*disturbance,2*control,2*control],'Curvature',[1,1],'EdgeColor',[0 0 0],'FaceColor','y');

axis([xi xf yi yf])

% axis square;
% axis equal tight
axis equal

set(gca,'fontsize',18)

% xlabel('$x$','Interpreter','latex','fontsize',30);

% ylabel('$\dot{x}$','Interpreter','latex','fontsize',30);

box on;

% if clusterFlag,
    saveas(h,['figure_',num2str(i),'_',num2str(control),'.fig'],'fig')	%Save the figure    
% end

if (i>1)
close(formerh);
end;

formerh=h;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PRINT ZONE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (current_num_points==previous_num_points || current_num_points==1)    %This is the stop condition. If the number of points in grid_points_base after the current cut
    break;                                                               %is equal to the number of points in the previous cut, then grid_points_base contains a
end;                                                                     %Safe Set.

previous_num_points=current_num_points;

runTimeSafeSet(i) = toc;

end;

% figure
% 
% plot(mem_points)
% 
% title('Here we show the number of points per iteration')
% 
% xlabel('Iteration');
% 
% ylabel('Number of points');
% 

% save(['safe_set_',num2str(epsilonC),'.mat'],'grid_points_base')
% save(['safe_set_sculpting_',num2str(epsilonC),'.mat'])

save(['safe_set_',num2str(e),'_',num2str(R),'_',num2str(control),'.mat'],'grid_points_base')
save(['safe_set_sculpting_',num2str(e),'_',num2str(R),'_',num2str(control),'.mat'])

end
