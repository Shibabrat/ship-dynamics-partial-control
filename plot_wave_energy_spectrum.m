%SCRIPT to plot the wave energy spectrum generating the rough seas for
%diffeent significant wave heights
%
% Shibabrat Naik
% Ref: http://blogs.mathworks.com/loren/2007/12/11/making-pretty-graphs/

% set(0,'defaultFigurePaperPositionMode', 'auto', ...
%     'DefaultAxesFontSize',20)

global H chi
chi=90*(pi/180);
figure('Units', 'pixels', ...
    'Position', [100 100 500 375]);
hold on;
box on;
for H=2:1:10

    %Generate rough sea waves 
    [omegaI,waveEnergySpect,~,~] = func_get_random_waves(0.1); %0.1 is arbitrary here
    hSpect = plot(omegaI,waveEnergySpect,'-r')
    
    hXLabel = xlabel('$\omega$ (rad/sec)','interpreter','latex')
    hYLabel = ylabel('$S$ (m$^2$-sec)','interpreter','latex')


end
hText  = text(1, 10, 'H_s = 2:1:10 m','FontSize', 14);
set(hSpect, 'LineWidth', 2);
set(gcf, 'PaperPositionMode', 'auto');
set([hXLabel, hYLabel], 'FontSize', 18);
set(gca,'Fontsize',14)
% print -depsc2 finalPlot1.eps
% close;

