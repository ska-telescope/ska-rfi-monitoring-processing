# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 11:08:17 2019

@author: f.divruno
"""

import numpy as np
import matplotlib.pyplot as plt

def smooth(x,window_len=11,window='hanning'):

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
#    return y
    return y[(int(window_len/2)-1):-int(window_len/2)] 

def find_RFI_mask_points(D,freq,k=5,test=0):
    '''
    Author: Federico Di Vruno
    Date: 10/10/2019
    
    This function looks for the frequency points in a datacube (time,freq) that do not have RFI
    components, this is accomplished using the standard deviation of the dataset 
    in the time axis.
    This frequency points can be later used to generate a mask of the "baseline"
    data interpolating with this points.
    IN:
        D: original data
        Freq: freq vector
        k: multiplier to threshold the data
        test: flag to run the built in test of the function.
    
    OUT:
        Frequency vector of the points without RFI detection.
        
    TODO: implement more complex statistics to detect RFI and not astronomical signals.
    '''  
    flags = np.zeros(np.shape(D))
    for i in range(len(D)):
        Dnorm = abs(D[i] - smooth(D[i],50))
        flags[i,:] = np.array(Dnorm >=0.6).astype('int')
     
    sum_flags = np.sum(flags,0)
    ind = np.where(sum_flags < len(flags)*0.001)[0]
    

    
    return freq[ind],flags


def occupancy3(freq,maskFreq,D,test =0):
    '''
        Author: Federico Di Vruno
        Date: 10/10/2019

        Calculates the percentage of the time that the power in a certain frequency is
        above a threshold.
        IN:
            freq:
            maskFreqs: array with frequencies without RFI as closely spaced as possible.
            D: data cube in linear units
            sigma_mult:
            test: flag to generate figures to see what the function is doing (BIT)
                
        OUT:
            Occupancy vector in % per frequency
            
        TODO:
            
    '''
    
    nTime = np.size(D,0)
    nFreq = len(freq)
    occupancy = np.zeros(nFreq)
    
    #for each timestep the flags are calculated and added to the occupancy vector.
    for i in range(nTime):
        D1 = D[i,:]
   
        mask_ave = np.interp(maskFreq,freq,D1) # resample the average vector in the mask points.
        mask_ave = np.interp(freq,maskFreq,mask_ave) #interpolate the mask back to the original sample number
        D2 = D1 - mask_ave # leave only the RFI

        Occup = np.array(D2 > 0).astype('int')         
        occupancy += Occup
        print(str(i) + ' of ' + str(nTime))
    
    if test:
        # in case needed the acerage plot and the mask can be plotted.
        # also plots the occupancy as final result.
        plt.figure(figsize = [15,10])
        ax = plt.axes()
        ax.plot(freq,10*np.log10(np.transpose(D1))) #original data
        ax.plot(freq,10*np.log10(np.transpose(D2))) #data without average
        ax.plot(freq,10*np.log10(mask_ave),label = 'mask std')
        ax.plot(freq,10*np.log10(D1),label = 'average')
        ax.scatter(maskFreq,10*np.log10(np.interp(maskFreq,freq,D1)),label = 'mask points in average')
        plt.xlim([freq[0],freq[-1]])
        plt.legend()
        plt.xlabel('Freq Mhz')
        plt.ylabel('log ADC units')
#        plt.savefig(outdir+ 'Average_occup_mask_'+ time_freq , dpi=600, bbox_inches='tight')

#        plt.figure(figsize = [15,10])
#        ax = plt.axes()
#        ax.plot(freq,10*np.log10(np.transpose(D1)))
#        plt.xlim([freq[0],freq[-1]])
#        plt.legend()
#        plt.xlabel('Freq Mhz')
#        plt.ylabel('log ADC units')
        
        
#        plt.figure(figsize = [15,10])
#        ax = plt.axes()
#        ax.plot(freq,occupancy*100/nTime)
#        plt.title('Occupancy')
#        plt.xlabel('Freq Mhz')
#        plt.ylabel('Occupancy % of time')
#        plt.savefig(outdir+ 'Occupancy_'+ time_freq , dpi=600, bbox_inches='tight')
        
    return occupancy*100/nTime



if __name__ == '__main__':
    #do integrated test, need a dataset
    import rfilib
    data = data_lin
    freq = freqs1
    freqRFI, flags = find_RFI_mask_points(data,freq,test=1)
    timeVect = np.linspace(0,len(data)-1,len(data))
    rfiLib.plot_spectrogram(timeVect,freq,data,'MeerKAT data')
    
    plt.figure()
    plt.plot(freq,ocup)
    
    plt.figure()
    plt.plot(freq, 10*np.log10(np.transpose(data[0:100,:])))
    plt.scatter(freqRFI, 10*np.log10(data_base))
    
    