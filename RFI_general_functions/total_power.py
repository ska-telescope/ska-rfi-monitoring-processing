# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 12:16:41 2019

@author: f.divruno
"""
import numpy as np
import matplotlib.pyplot as plt

def total_power(time,freq,data,fo,B,Tlow=0,Thigh=1e100,plot_flag=0):
    #Calculates the total power in the specified frequency band, can include a threshold calculation and 
    #a flag  to plot the data.
    # input data should be in V^2, or power.
    
    fmin = fo-B/2
    fmax = fo+B/2
    
    fstep = freq[1] - freq[0]
    
    ind1 = int((fmin-freq[0])/fstep)
    ind2 = int((fmax-freq[0])/fstep)
    
    total_power = np.sum(data[:,ind1:ind2]**2,1)
    
    # include some threshold limitations
    if (Tlow!=0 and Thigh!=1e100):
        mask = np.ones(len(total_power), dtype=bool)
        for i in range(np.size(total_power)):
            if (total_power[i] < Tlow or total_power[i] > Thigh):
                mask[i] = False
        Data2 = total_power[mask]
        time = time[mask]
        # include some threshold limitations
    
    
    if plot_flag ==1:
        plt.figure()
        plt.plot(10*np.log10((Data2)))
        plt.title('freq = ' + str(fo) + ' MHz, B = ' + str(B) + ' MHz')
        plt.grid(True,'both')
    return time,Data2
