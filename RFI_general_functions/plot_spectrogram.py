# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 12:20:55 2019

@author: f.divruno
"""
import numpy as np
import matplotlib.pyplot as plt

def plot_spectrogram(time,freq,data,title,Fstart=0,Fstop=0,Tstart=0,Tstop=0):
    if Fstop==0:
        Fstop = freq[-1]
        
    if Tstop==0:
        Tstop = time[-1]
    
    # mask the frequency
    mask1= np.array(freq>Fstart)
    mask2 = np.array(freq<Fstop)
    mask3 = mask1 & mask2
    freq = freq[mask3]
    data = data[:,mask3]   
    
    #mask the time
    mask1= np.array(time>Tstart)
    mask2 = np.array(time<Tstop)
    mask3 = mask1 & mask2
    time = time[mask3]
    data = data[mask3,:]   
    
    
    fig = plt.figure(figsize=(20,12))
    ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
    
    left = freq[0]
    right = freq[-1]
    bottom = time[0]
    top = time[-1]
    data_log = (10*np.log10(data))
    cim = ax.imshow( data_log, origin='lower', interpolation='nearest', cmap= 'jet', extent=(left, right, bottom, top), )
     
    ax.set_aspect(abs(right-left) / abs(top-bottom))
    plt.xlabel('MHz')
    plt.ylabel('time')
    plt.title(title)   

