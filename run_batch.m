%	Script to run batch mode for asymptotic safe set 
% 	Authors: 
%			Shibabrat Naik, Shane Ross
%	Last modified:
%           10 May 2015
%           03 May 2015
%			01 April 2015


%% Roll model for ship dynamics
% global epsilonC
% for epsilonC = 1.95:0.05:3
%     fileName = ['diary_',num2str(epsilonC)];
%     diary(fileName)
%     safe_set_finder(1001,0)
%     diary off
%     close all
% end


%% Roll-pitch model for ship dynamics
% [~,grid_points_base] = func_get_initial_set_rect(0.253,1.6,1001,1001);
% func_get_image_rect(0.253,1.6,grid_points_base)
% safe_set_finder(500,0)
global control disturbance
disturbance = 0.1;
for control = 0.0:0.0025:disturbance
    fileName = ['diary_',num2str(control),'.txt'];
    diary(fileName)
    safe_set_finder(500,0)
    diary off
    close all
end