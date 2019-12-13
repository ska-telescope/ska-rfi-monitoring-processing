# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 22:58:52 2019
    Flagging function
    
@author: f.divruno
"""
import numpy as np
import matplotlib.pyplot as plt

def flagging(D, thr = 0.1, time=[], freq=[], test = 0):
    '''
        Simple threshold flagging.
        input:
            D in log units
            Thr in log units
            time and freq are for plotting in case of test=1
        output:
            flags: array of same size as D with 1 where a flag whas detected
            occupancy: sum of all the flags in the time vector.
    '''
    def smooth(x,window_len=20):
        s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
        w=np.hanning(window_len)
        y = np.convolve(w/w.sum(),s,mode='valid')
        return y[(int(window_len/2)-1):-int(window_len/2)] 

    flags = np.zeros(np.shape(D))
    for i in range(len(D)):
        Dnorm = abs(D[i] - smooth(D[i],50))
        flags[i,:] = np.array(Dnorm >=thr).astype('int')
     
    occupancy = np.sum(flags,0)/len(flags)*100
    
    if test :
        from rfiLib.plot_spectrogram import plot_spectrogram
        
        plt.figure()
        plt.plot(Dnorm)
        plt.plot([0,len(Dnorm)],[thr,thr])
        
        if time==[]:
            time = np.arange(0,len(D))
        if freq==[]:
            freq = np.arange(0,np.size(D,1))
            
        plot_spectrogram(time,freq,10**(D/10),'Original data')   
        plot_spectrogram(time,freq,flags,'flags')   
        
        
    return flags, occupancy



