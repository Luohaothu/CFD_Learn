# -*- coding: utf-8 -*-

"""
CFD Solver of subsonic-supersonic isentropic nozzle flow

Method : Maccormack's Technique
"""

import numpy as np
from pylab import *

n = 61 #node numbers
L = 3.0  #total length
rou = np.zeros(n)       #density
v   = np.zeros(n)       #velocity
T   = np.zeros(n)       #tempeture
A   = np.zeros(n)       #section area
x   = np.linspace(0, L, n)  #node ordinate
dx  = L / (n - 1)       #space steplength
dt_l= np.zeros(n)       #time steplength list
C   = 0.5               #Courant number
STEP_N = 1000           #step number
t   = 0                 #total time
gama= 1.4               #gama

A = 1 + 2.2*(x - 1.5)**2
rou = 1 - 0.314 * x
T = 1 - 0.2314 * x
v = (0.1 + 1.09 * x) * np.sqrt(T)

##plot init value
#subplot(2,2,1)
#plot(x, A)
#
#subplot(2,2,2)
#plot(x, rou)
#
#subplot(2,2,3)
#plot(x, T)
#
#subplot(2,2,4)
#plot(x, v)
#
#show()

rou_hist = []
v_hist = []
T_hist = []
M_hist = []
dvdt_hist = []
drdt_hist = []
dTdt_hist = []
#Calculation
for i in xrange(STEP_N):
    rou_hist.append(rou[n/2])
    v_hist.append(v[n/2])
    T_hist.append(T[n/2])
    M_hist.append(v[n/2] / sqrt(T[n/2]))
    dvdt_hist.append(log10(abs(np.max(dvdt))))
    drdt_hist.append(log10(abs(np.max(drdt))))
    dTdt_hist.append(log10(abs(np.max(dTdt))))
    dt_l = C * dx / (abs(v) + np.sqrt(T))
    dt = np.min(dt_l)
    t += dt
    
    #predict-step
    drdt1 = -v * (roll(rou, -1) - rou) / dx - rou * (roll(v, -1) - v) / dx \
           -rou * v * (np.log(roll(A, -1)) - np.log(A)) / dx
    dvdt1 = -v * (roll(v, -1) - v) / dx - 1 / gama * ((roll(T, -1) - T) / dx \
           + T / rou * (roll(rou, -1) - rou) / dx)
    dTdt1 = -v * (roll(T, -1) - T) / dx - (gama - 1) * T *((roll(v, -1) - v)/dx\
           + v * (np.log(roll(A, -1)) - np.log(A)) / dx)
    
    drdt1[-1] = 0
    dvdt1[-1] = 0
    dTdt1[-1] = 0
    drdt1[0] = 0
    dvdt1[0] = 0
    dTdt1[0] = 0

    
    rou_pre = rou + drdt1 * dt
    v_pre   = v + dvdt1 * dt
    T_pre   = T + dTdt1 *dt
        
    #correct-step
    drdt2 = -v_pre * (rou_pre - roll(rou_pre, 1)) / dx - rou_pre * (v_pre - roll(v_pre, 1)) / dx \
            -rou_pre * v_pre * (np.log(A) - np.log(roll(A, 1))) / dx
    dvdt2 = -v_pre * (v_pre - roll(v_pre, 1)) / dx - 1 / gama * ((T_pre - roll(T_pre, 1)) / dx \
            +T_pre / rou_pre * (rou_pre - roll(rou_pre, 1)) / dx)
    dTdt2 = -v_pre * (T_pre - roll(T_pre, 1)) / dx - (gama - 1) * T_pre * ((v_pre - roll(v_pre, 1)) / dx \
            +v_pre * (np.log(A) - np.log(roll(A, 1))) / dx)
    
    drdt2[-1] = 0
    dvdt2[-1] = 0
    dTdt2[-1] = 0
    drdt2[0] = 0
    dvdt2[0] = 0
    dTdt2[0] = 0


    drdt = (drdt1 + drdt2) / 2
    dvdt = (dvdt1 + dvdt2) / 2
    dTdt = (dTdt1 + dTdt2) / 2
    
    rou = rou + drdt * dt
    v   = v   + dvdt * dt
    T   = T   + dTdt * dt
    
    v[0] = 2 * v[1] - v[2]
    v[-1] = 2 * v[-2] - v[-3]
    rou[0] = 1.0
    rou[-1] = 2 * rou[-2] - rou[-3]
    T[0] = 1.0
    T[-1] = 2 * T[-2] - T[-3]
    
    if i / 10 == 0 :
        print "\n    t \t   dt \t  drdt\t  dvdt\t  dTdt\t   rou\t    v \t    T \n"
    print "%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t\n" %(t, dt, np.average(drdt), np.average(dvdt), \
            np.average(dTdt), np.average(rou), np.average(v), np.average(T))
    
subplot(2,3,1)
plot(x, v/np.sqrt(T))
title('Mach number')
subplot(2,3,2)
plot(np.array(range(STEP_N)), np.array(rou_hist))
title('rou history data')
subplot(2,3,3)
plot(np.array(range(STEP_N)), np.array(v_hist))
title('velocity history data')
subplot(2,3,4)
plot(np.array(range(STEP_N)), np.array(T_hist))
title('temperture history data')
subplot(2,3,5)
plot(np.array(range(STEP_N)), np.array(M_hist))
title('Mach number history data')
subplot(2,3,6)
plot(np.array(range(STEP_N)), np.array(drdt_hist))
plot(np.array(range(STEP_N)), np.array(dvdt_hist))
plot(np.array(range(STEP_N)), np.array(dTdt_hist))
title('partial derivative history data')
show()