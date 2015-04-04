%	Script to run batch mode for asymptotic safe set 
% 	Authors: 
%			Shibabrat Naik, Shane Ross
%	Last modified:
%			01 April 2015
N = 1001;

global epsilonC
for epsilonC = 1.9:0.05:3
    fileName = ['diary_',num2str(epsilonC)];
    diary(fileName)
    asymptotic_sculpting_finder(N)
    diary off
    close all
end

