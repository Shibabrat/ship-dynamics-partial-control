function h = plot_nondim_potential()
        
    set(0,'defaulttextinterpreter','latex')
    set(0,'DefaultAxesFontName','Times New Roman');
    set(0,'DefaultTextFontName','Times New Roman');
    set(0,'DefaultAxesFontSize',26);
    set(0,'DefaultTextFontSize',20);
    
    h = ezplot(@(x)nondim_potential(x),[-1 1]);
    title([])
    xlabel('$\phi$')
    ylabel('$V(\phi)$')
    box on
    grid on
    
    line(0.88*ones(1000,1),linspace(-0.3,3,1000),'Color','m')
    line(-0.88*ones(1000,1),linspace(-0.3,3,1000),'Color','m')
    text(-0.85,-0.01, ...
        ' $\longleftarrow$ Roll angle of vanishing stability $\longrightarrow$ ')

end
function V = nondim_potential(x)

    %Parameters for the Edith Terkol
    omegaN = 0.62;
    c1Bar = omegaN^2;
    c2Bar = 0.1296/c1Bar;
    c3Bar = 1.0368/c1Bar; 
    c4Bar = -4.059/c1Bar;
    c5Bar = 2.4052/c1Bar;
    
    V = 0.5*(x.^2) + (c2Bar/3)*(abs(x).*(x.^2)) + ...
        (c3Bar/4)*(x.^4) + (c4Bar/5)*(abs(x).*(x.^4)) + ...
        (c5Bar/6)*(x.^6);

end