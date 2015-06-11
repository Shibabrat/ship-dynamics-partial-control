# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 01:44:30 2015

Script to test the event option for ODE integration using PYCSE module

@author: shibabrat
"""

#import warnings
import math 
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import sys
sys.path.append('/Library/Python/2.7/site-packages/')
from pycse import *


#Setting the parameters for the model
omegaN = 0.62
omegaE = 0.527
omegaBar = omegaE/omegaN
b1 = 0.0043/omegaN
b2 = 0.0225
c2=0.1296/(omegaN**2)
c3=1.0368/(omegaN**2)
c4 = -4.059/(omegaN**2)
c5 = 2.4052/(omegaN**2)

#Forcing due to regular seas
alpha0 = 0.73;
H = 9.84;
lam = 221.94;
HBar = alpha0*math.pi*(H/lam);

def func_odes_boat_roll(x,t):
    
    phi, pPhi = x   
    dPhidt = pPhi
    dPphidt = - phi - c2*np.abs(phi)*phi - c3*np.power(phi,3) -\
        c4*(np.abs(phi)*np.power(phi,3)) - c5*np.power(phi,5) -\
        b1*pPhi - b2*(np.abs(pPhi)*pPhi) + HBar*np.sin(omegaBar*t)

#    print phi
#    print pPhi

    return [dPhidt, dPphidt]

def func_event_capsize(x, t):
    isterminal = False
    direction = 1
    value = np.abs(x[0]) - 0.88
    return value, isterminal, direction
    
def func_roll_model():
    x0 = 0.15
    y0 = 0.57
    #x_init = np.array([x0,y0])
    x_init = [x0,y0]
    time = np.linspace(0.0,4*((2*np.pi)/omegaBar),1000)
    #x_sol = odeint(func_odes_boat_roll,x_init,time) 
    
    t_sol, x_sol, te, xe, ie = odelay(func_odes_boat_roll,x_init,time, \
                                  events=[func_event_capsize])

    #plt.plot(time,x_sol[:,0],'-b')
    #plt.hold(True)
    #plt.plot(time,0.88*np.ones(1000),'--r',time,-0.88*np.ones(1000),'--r')

    plt.plot(t_sol,x_sol[:,0],'-b')
    plt.hold(True)
    plt.plot(time,0.88*np.ones(1000),'--r',time,-0.88*np.ones(1000),'--r')
    plt.plot(te,xe[:,0],'om')    
    plt.show()
    
func_roll_model()









