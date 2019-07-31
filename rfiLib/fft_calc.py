# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 12:21:44 2019

@author: f.divruno
"""
import numpy as np
import matplotlib.pyplot as plt

def fft_calc(signal, fs, power_calc=1 ,plot_figs=0):
    # if power_calc = 1: Returns the power spectrum in real frequencies and the freq vector in MHz.
    # if power_calc = 0: Returns the complex voltage spectrum and the freq vector in MHz.
    
    max_length =  10*131072
    try:
        N = np.size(signal,1)
        N_files = np.size(signal,0)
    except:
        N = len(signal)
        N_files = 1        

    
    
    freq = (np.fft.fftfreq(N, d=1/fs))
    if power_calc ==1:
        data_fft = partition(signal,max_length,np.fft.fft,power_spectrum)/N
        freq = freq[0:int(N/2)]
    else:
        data_fft = partition(signal,max_length,np.fft.fft)/N
    
    plot_all = 0 # this is for debugging, if set to 1 generates 1 plot for each file.
    if plot_all ==1:
#        N_files = 1 #for  debugging
        plt.figure()
        for i in range(N_files):
            D = abs(data_fft)**2
            plt.plot(freq/1e6,10*np.log10(D))    
        plt.title('power spectrum' )

    if plot_figs==1:
        if N_files >1:
            ave = np.average(abs(data_fft)**2,0)
        else:
            ave = np.transpose(abs(data_fft)**2)
        plt.figure()
        plt.plot(freq/1e6,10*np.log10(ave))    
        plt.title('Average of all the captures')
    return [freq/1e6,data_fft]

def idle_fun(val):
    return val

def partition(signal,max_len, func1, func2= idle_fun):
    # This function breaks a very long dataset in chunks to process func1 and func2 (more could be added)
    # This is used to calculate for example the FFT of a m x n samples, the m dimension is 
    # divided in "max_len/n" chunks and the result stiched in the end.

    try:
        m = np.size(signal,0)
        n = np.size(signal,1)
    except:
        m = 1
        n = len(signal)
    data = np.ndarray([m,int(n/2)])
    if (m*n >  max_len) & (m>1): 
        steps = int(m*n/max_len)
        N_step = int(m/steps)
        for i in range(steps):
            data[int(i*N_step):int((i+1)*N_step),:] = func2(func1(signal[int(i*N_step):int((i+1)*N_step),:]))
            print(str(i+1)+' of '+str(steps)+' steps in partition')
        if m > (i+1)*N_step:
            data[int((i+1)*N_step)::,:] = func2(func1(signal[int((i+1)*N_step)::,:]))
    else:
        data= func2(func1(signal))
    return data


def power_spectrum(data):
    try:
        N = np.size(data,1) 
        P = np.abs(data[:,0:int(N/2)]).astype('float32')**2
    except:
        N = len(data)
        P = np.abs(data[0:int(N/2)]).astype('float32')**2
    
    
    return P