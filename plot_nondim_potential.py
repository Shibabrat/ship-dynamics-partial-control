# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 08:59:19 2015

@author: shibabrat
"""

import numpy as np
import matplotlib.pyplot as plt

#%% Non-dimensional potential for roll model in ship dynamics
def nondim_potential():
    phi = np.linspace(-1.5,1.5,1000,endpoint=True)

    omegaN = 0.62
    c1Bar = np.power(omegaN,2.0)
    c2Bar = 0.1296/c1Bar
    c3Bar = 1.0368/c1Bar 
    c4Bar = -4.059/c1Bar
    c5Bar = 2.4052/c1Bar
    nondimPot = 0.5*np.power(phi,2.0) + \
        (c2Bar/3.0)*(abs(phi)*np.power(phi,2.0)) + \
        (c3Bar/4.0)*np.power(phi,4.0) + \
        (c4Bar/5.0)*(abs(phi)*np.power(phi,4.0)) + \
        (c5Bar/6.0)*np.power(phi,6.0)        
    
    plt.plot(phi,nondimPot,'--b')
    plt.grid()   
    plt.xlabel('$\phi$',fontsize = 30)    
    plt.ylabel('$V(\phi)$',fontsize = 30)    
    
    plt.show()
    
nondim_potential()


   
        
    
    
