function asymptotic_sculpting_finder(N,clusterFlag)
% This algorithm computes the Asymptotic Sculpting Safe Set from a known safe set 
% stored in the file safe_set. It also needs the forward iteration
% of the points that form the safe set, which are stored in the file 
% forward_iteration. To make it work properly it has to be given the 
% following parameters: 
% control: the value of the bound of the control.
% disturbance: the bound of the value of the disturbance.
% xi,yi,xf,xf: the corners of the square where we want to find the safe set.
% N: the resolution of the grid, being the total number of points N*N.
%
% It uses three different grid if points:
%
% -The array grid_points_asymptotic safe_set that will contain the initial
% set from which we will start sculpting and the succesive steps until the algorithm 
% converges to the desired Asymptotic safe set.
% -The array grid_points_iteration that will contain the forward iteration
% of all the points in the set grid_points_asymptotic.
% -The array grid_points_reachable that will contain all the points that can
% be reach from the set grid_points_iteration for any possible combination 
% of disturbance plus control.
%
% Adapted and modified by Shibabrat Naik, Shane Ross
% Date: 06 March 2015

% clear all
% clc
% close all

format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CONFIGURATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% control=0.0475;    %Here we specify the bound of the control that we are going to use.

% disturbance=0.08;    %Here we specify the bound of the disturbance that we are going to use.


% xi=-2;     %Here are the corners of the square where we are going to
% yi=-2;     %compute the Asymptotic Safe Set.

% xf=2;
% yf=2;

global epsilonC
% epsilonC = 2.2;
epsilonD = 3;
H = 4.94;
HBar = 0.73*pi*(H/221.94);
control = epsilonC*HBar
disturbance = epsilonD*HBar
safeRatio = control/disturbance
omegaN = 0.62;

xi=-0.88;    %Here are the corners of the box where we are going to 
yi=-0.52;    %compute the Safe Sets and  the Asymptotic Safe Set.

xf=0.88;
yf=0.52;

if nargin == 0 	%Default case
	N=3001;    	%Here we specify the resolution of the grid, for example if we choose
           		%N=6001, we will have 6001x6001=36012001 points.
end

save(['asymptotic_sculpting_params_',num2str(epsilonC),'.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CONFIGURATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%====LCC Cluster specific===============%
if clusterFlag,
    load(['~/partial-control-roll-capsize/caseD/Fig5c-caseG-control-disturbance-I/safe_set_',num2str(epsilonC),'.mat']);
    load('~/partial-control-roll-capsize/caseD/initial_set_image.mat');
else
    load(['safe_set_',num2str(epsilonC),'.mat'])
    load initial_set_image
end
%=======================================%

[indexes_reachable]=index(N,grid_points_base,disturbance,control);    %Here we compute the relative indexes needed to carry out the fatten operation.    

previous_num_points=0;

fig_points=zeros(N*N,2);    %In this array we will store the points that remain after each cut in order to plot them.

grid_points_asymptotic=grid_points_base;


for i=1:100
    
fprintf('\n\nThis is the cut number %d\n\n',i);
    
disp('Step 1');
grid_points_iteration=iteration(N,grid_points_asymptotic,grid_points_image,xi,xf,yi,yf);    %Here we compute the forward iteration of the set that remains after every cut and is stored in grid_points_asymptotic.    
disp('Step 2');   
boundary_index=iteration_boundary(N,grid_points_iteration);    %Here we compute the boundary of the set grid_points_iteration.
disp('Step 3');
grid_points_reachable=iteration_fatten(N,grid_points_iteration,boundary_index,indexes_reachable);    %Here we compute the fatten operation of the set grid_points_iteration.
disp('Step 4');
grid_points_asymptotic=intersection(N,grid_points_base,grid_points_reachable);    %Here we compute the intersection between the sets grid_points_reachable and the original safe set that we have in grid_points_base.
disp('Step 5');



current_num_points=1;    %We use this variable to count the numbers of points that remain after every cut.

fig_points=zeros(N*N,2);    %In this array we will store the points that remain after each cut in order to plot them.

for j=1:N    %Here we count the number of points that remain after every cut. 
     for k=1:N
             if(grid_points_asymptotic(1,3,j,k)==1)
                 fig_points(current_num_points,:)=grid_points_asymptotic(1,1:2,j,k);
                 current_num_points=current_num_points+1;
             end;
     end;
end;

mem_points(i)=current_num_points    %In this array are stored the succesive number of points that remain in each iteration.

fig_points=fig_points(1:current_num_points-1,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PRINT ZONE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=figure;

if clusterFlag,
    set(h,'Visible','off')		%When display is not allowed in compute clusters
end

plot(fig_points(:,1),fig_points(:,2),'r.')
% plot(fig_points(:,1),fig_points(:,2),'.b')


iteration_number=num2str(i);

title(['Iteration number ' iteration_number ' of the Asymptotic Sculpting Algorithm'])


rectangle('Position',[xi+4*disturbance,yi+8*disturbance,2*disturbance,2*disturbance],'Curvature',[1,1],'EdgeColor',[0 0 0],'FaceColor','g');
rectangle('Position',[xi+4*disturbance+disturbance-control,yi+4*disturbance,2*control,2*control],'Curvature',[1,1],'EdgeColor',[0 0 0],'FaceColor','y');

% axis square;
% axis equal tight
axis equal

axis([xi xf yi yf]);

xlabel('$x$','Interpreter','latex');

ylabel('$\dot{x}$','Interpreter','latex');

box on;

if clusterFlag,
    saveas(h,['figure_',num2str(i),'.fig'],'fig')	%Save the figure    
end

if (i>1)
close(formerh);
end;

formerh=h;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PRINT ZONE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if current_num_points==previous_num_points    %This is the stop condition. If the number of points in grid_points_asymptotic after the current cut
    break;                                    %is equal to the number of points in the previous cut, then grid_points_asymptotic contains an 
end;                                          %Asymptotic Safe Set.

previous_num_points=current_num_points;

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


save(['asymptotic_safe_set_',num2str(epsilonC),'.mat'],'grid_points_asymptotic')
save(['asymptotic_sculpting_',num2str(epsilonC),'.mat'])


end

