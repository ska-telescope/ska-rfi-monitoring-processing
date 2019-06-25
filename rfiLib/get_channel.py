# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 12:16:41 2019

@author: f.divruno
"""
import numpy as np
import matplotlib.pyplot as plt

def get_channel(freq,data,fo,B,plot_flag=0):
    #Calculates the total power in the specified frequency band, can include a threshold calculation and 
    #a flag  to plot the data.
    # input data should be in V^2, or power.
    
    fmin = fo-B/2
    fmax = fo+B/2
    
    fstep = freq[1] - freq[0]
    
    N = np.size(data[1,:])
    M = np.size(data[:,1])
    fs = (freq[1] - freq[0])*N*2
    
    ind1 = int((fmin-freq[0])/fstep)
    if ind1<1:
        ind1=0
        
    ind2 = int((fmax-freq[0])/fstep)
    if ind2>len(freq):
        ind2 = len(freq)
    
    fd_chan = np.concatenate(np.zeros([M,ind1-1]), data[:,ind1:ind2] , np.zeros([M,ind2-1]))
    fd_chan = np.concatenate(fd_chan[N::-1],fd_chan[0:N])
       
    td_chan = np.real(np.fft.ifft(fd_chan))
    time = np.linspace(0,1/fs*(N-1),N)    
    
    if plot_flag ==1:
        plt.figure()
        plt.plot(10*np.log10((fd_chan)))
        plt.title('freq = ' + str(fo) + ' MHz, B = ' + str(B) + ' MHz')
        plt.grid(True,'both')

        plt.figure()
        plt.plot(time,td_chan)
        plt.title('Time signal in ' + str(fo) + ' MHz, B = ' + str(B) + ' MHz')
        plt.grid(True,'both')

    
    return time,td_chan
