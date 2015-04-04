%	Script to run batch mode for asymptotic safe set 
% 	Authors: 
%			Shibabrat Naik, Shane Ross
%	Last modified:
%			01 April 2015

global epsilonC
for epsilonC = 1.95:0.05:3
    fileName = ['diary_',num2str(epsilonC)];
    diary(fileName)
    safe_set_finder(1001,0)
    diary off
    close all
end