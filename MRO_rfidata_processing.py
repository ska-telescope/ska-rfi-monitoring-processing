# -*- coding: utf-8 -*-
"""
Created on Wed May 22 02:42:48 2019

@author: F.Divruno
"""
import matplotlib.pyplot as plt
import matplotlib
import rfiLib as RFI
import numpy as np
import sys
import os, os.path
import subprocess

#Commands to run in Ubuntu server:

# tmux to use a deaachable console.
# python3 to load the interactive python (so that the data remains in memory)
# subprocess.run(["git","pull"])
# exec(open("./MRO_rfidata_processing.py").read())
# subprocess.run(["git","pull"]) to use githul pull
# exec(open("./MRO_rfidata_processing.py").read()) to execute the script from the interactive python.



maskFreq =[20, 22, 25, 27.5, 31, 34.1, 37.2, 47.2,
           55.8, 59.5, 60, 70, 75.6, 79, 80, 83, 87,
           92, 94, 103, 110, 117, 119.5, 133, 136,
           140, 144, 149, 151, 155, 159,
           159.4, 162, 163, 165, 167, 169, 171, 173,
           175, 177, 179, 181, 183, 185, 187,189, 191, 
           193, 195, 197, 199, 201, 203, 205, 207, 209,
           211, 213, 215, 217, 219, 221, 223, 225, 227,
           229, 231, 233, 235, 237, 239,240, 242,246,
           253, 259, 264, 269.2, 271, 272, 274, 276,
           278, 280, 282, 284, 286, 288, 290, 292, 294,
           296, 298, 300, 302, 304, 306, 308, 310, 312,
           314, 316, 318, 320, 322, 324, 326, 328, 330,
           332, 334, 336, 338, 340, 342, 344, 346, 348,
           350, 352, 354, 356, 358, 365, 370, 375, 380,
           382, 384, 386, 388, 390, 392, 394, 396, 398,
           400]

def occupancy(freq,maskFreq,D,sigma_mult,plot_figs =0):
    '''
        Calculates the percentage of the time that the power in a certain frequency is
        above a threshold.
        inputs:
                freq:
                maskFreqs: array with frequencies without RFI as closely spaced as possible.
                D:
                sigma_mult:
                plot_figs:
    
    '''
    
    nTime = np.size(D,0)
    nFreq = len(freq)
    ave = np.average(D,axis=0)
    sigma = np.std(D,axis=0)
    occupancy = np.zeros(nFreq)
    mask_ave = np.interp(maskFreq,freq,ave)
    mask_ave = np.interp(freq,maskFreq,mask_ave)
    mask_sigma = np.interp(maskFreq,freq,sigma)
    mask_sigma = np.interp(freq,maskFreq,mask_sigma)
    
    
    
    mask = mask_ave+mask_sigma*sigma_mult
    
    #for each timestep the flags are calculated and added to the occupancy vector.
    navg = 10
    nTime = int(nTime/navg) # averages every 10 time samples
    for i in range(nTime):
        D_avg = np.average(D[i*navg:(i+1)*navg,:],0)
        occupancy += np.array(D_avg > mask).astype('int')
        
    # in case needed the acerage plot and the mask can be plotted.
    # also plots the occupancy as final result.
    if plot_figs:
        plt.figure(figsize = [15,10])
        ax = plt.axes()
        ax.plot(freq,10*np.log10(mask_sigma),label = 'mask')
        ax.plot(freq,10*np.log10(mask_sigma*sigma_mult),label = 'mask * multiplier')
        plt.legend()
        plt.xlabel('Freq Mhz')
        plt.ylabel('log ADC units')
        
        
        plt.figure(figsize = [15,10])
        ax = plt.axes()
        ax.plot(freq,10*np.log10(mask),label = 'mask')
        ax.plot(freq,10*np.log10(ave),label = 'average')
        ax.scatter(maskFreq,10*np.log10(np.interp(maskFreq,freq,ave)),label = 'mask points in average')
        plt.xlim([freq[0],freq[-1]])
        plt.legend()
        plt.xlabel('Freq Mhz')
        plt.ylabel('log ADC units')
        plt.savefig(outdir+ 'Average_occup_mask_'+ time_freq , dpi=600, bbox_inches='tight')
        
        plt.figure(figsize = [15,10])
        ax = plt.axes()
        ax.plot(freq,occupancy*100/nTime)
        plt.title('Occupancy')
        plt.xlabel('Freq Mhz')
        plt.ylabel('Occupancy % of time')
        plt.savefig(outdir+ 'Occupancy_'+ time_freq , dpi=600, bbox_inches='tight')
        
    return occupancy*100/nTime






# Clear the screen
os.system('clear')

if os.name == 'posix':
    print('Using Ubuntu settings')
    ubuntu = 1
if os.name == 'nt':
    print('Using Windows settings')
    ubuntu = 0
    
#ubuntu = int(input('Select Ubuntu (1) or Windows (0) : '))

if ubuntu:
# in ubuntu
    matplotlib.use('Agg')
    plt.close('all')
    outdir = r'/mnt/data/MRO_rfidata/MRO_rfidata_19_05_12/output/'
    indir = r'/mnt/data/MRO_rfidata/MRO_rfidata_19_05_12'
else:
# in windows
    indir = 'C:\\Users\\F.Divruno\\Dropbox (SKA)\\14- RFI environment\\01- Australia\\rfidata_mro\\'
    outdir = 'C:\\Users\\F.Divruno\\Dropbox (SKA)\\14- RFI environment\\01- Australia\\rfidata_mro\\results\\'



#%% Read the files from disk
    
read_files = int(input('---Input data--\n\nData in memory (0)\nRead saved npz file (1)\nRead Fits files (2)\n Selection: '))
if read_files==2:
    print('reading files...')
    [freq,data] = RFI.read_MRO_data(indir,outdir)
    # Save the full data loaded in one file:
    if int(input('Save npz file?\nYes(1) No(0): ')):
        np.savez_compressed(outdir + r'MRO_rfidata_full', freq=freq, data=data)
    
elif read_files==1:
    try:
        print('Loading the complete set of saved data...')
        Aux = np.load(outdir+r'MRO_rfidata_full.npz')
        data = Aux['data'] #in V**2 originally, scaled to get to dBm, a guess on the calibration factor.
        freq = Aux['freq'] # in MHz   
    except:
        print('Failed...')
        print('Loading the saved data by files one by one...')
        files = os.listdir(outdir)
#        print(files)
        data = np.zeros([0,29801]).astype('float32')
        for i in range(len(files)): # comment for debug
#        for i in range(6):
            if os.path.splitext(files[i])[1] == '.npz':
                print(files[i])
                Aux = np.load(outdir+files[i])
                data = np.concatenate((data,Aux.get('data_file')),0) #in V**2 originally, scaled to get to dBm
                freq = Aux['freq'] # in MHz   
        np.savez_compressed(outdir + r'MRO_rfidata_full', freq=freq, data=data)
else:
    print('Data already in memory:\nData = %d lines x %d freq points' % (len(data),len(data[0])))
    
#%% TAke a subset of the frequency range
t_step = 1200/680
tmax = (len(data)-1)*t_step
time1 = np.linspace(0,tmax,len(data))
    
selection = input('\n\nSelect action:\
                    \n1: change frequency and time range\
                    \n2: Histogram\n3: Integrated Power\
                    \n4: Average\n5: Percentiles\
                    \n6: Occupancy\
                    \n7: Cumulative distribution function\
                    \n0: exit\
                    \n\nSelect: ')
init = 1    
while selection != '0':
    if init: 
        selection = '1'
        init = 0
    if selection == '1': #Change freq range and time
        fmin_select = input('Minimum freq to analyze (MHz) (enter for fmin): ')
        fmax_select = input('Maximum freq to analyze (MHz) (enter for fmax): ')
        timestart = input('Start time in sec (enter for 0): ')
        timestop = input('End time in sec (enter for tmax): ')

        if fmin_select == '':
            fmin =  freq[0]
        else:
            fmin = float(fmin_select)

        if fmax_select == '':
            fmax =  freq[-1]
        else:
            fmax = float(fmax_select)
               
        if timestart == '':
            tmin =  0
        else:
            tmin = int(timestart)
        if timestop == '':
            tmax =  time1[-1]
        else:
            tmax = int(timestop)         
        print('Slicing data')
        i_f_start = int(np.where(freq>=fmin)[0][0])
        i_f_stop = int(np.where(freq<=fmax)[0][-1])
        i_t_start = int(np.where(time1>=tmin)[0][0])
        i_t_stop = int(np.where(time1<=tmax)[0][-1])
        
        Daux = data[:,i_f_start:i_f_stop]/6.5-200 # scaling to match the levels calculated by gianni
        D = Daux[i_t_start:i_t_stop,:]
        print('data sliced')
        freqs = freq[i_f_start:i_f_stop]
        time = time1[i_t_start:i_t_stop]
        print('freq and time sliced')
        Pow = 10*np.log10(np.sum(10**(D/10),1))       
        time_freq = str(int((tmax-tmin)/3600))+'hs_'+str(int(fmin)) + 'to'+str(int(fmax))+'MHz'

        print('Done')

    if selection == '2': #Histogram
        print('Calculating histogram of the integrated power in %d to %d MHz' %(fmin,fmax))
        plt.figure(figsize=(30,20),dpi=600)
        plt.hist(Pow,500, linewidth=2)
        plt.xlabel('Ampl [dBm]')
        plt.grid()
        plt.savefig(outdir+ 'power_histogram_'+ time_freq , dpi=600, bbox_inches='tight')

    if selection == '3': #Total power
        #timestart = input('Start time in sec (enter for 0): ')
        #timestop = input('End time in sec (enter for tmax): ')

        print('Calculating total power...')
        
#        if timestart == '':
#            tmin =  0
#        else:
#            tmin = int(timestart)
#
#        if timestop == '':
#            tmax =  time[-1]
#        else:
#            tmax = int(timestop)
#        
#        time2 = time[(time>tmin) & (time<=tmax)]
#        Pow2 = Pow[(time>tmin) & (time<=tmax)]
            
        plt.figure(figsize=(30,20),dpi=600)
        plt.plot(time,Pow, linewidth=2)
        plt.xlabel('time [seconds]')
        plt.ylabel('Ampl [dBm]')
        title = 'total_power_'+ time_freq
        plt.title(title)
        plt.grid()
        plt.savefig(outdir+title , dpi=600, bbox_inches='tight')
        print('Done')


    if selection == '4': #Average
        print('Calculating Average')
        ave = np.average(D,0)            
        plt.figure(figsize=(15,10))
        plt.plot(freqs,ave, linewidth=2)
        plt.ylabel('Ampl [dBm]')
        plt.xlabel('freq [MHz]')
        plt.grid()
        title = 'Average_'+ time_freq
        plt.title(title)
        
        plt.savefig(outdir+ title, dpi=600, bbox_inches='tight')


    if selection == '5': #Percentiles00
        perc = input('Percentile = ')

        print('Calculating %s percentile' % (perc))
        
        title = 'MRO data '+ time_freq
        perc = int(perc)
        RFI.plot_percentile(freqs,D,perc,outdir,'dBm',title)

 
    if selection == '6': #Occupancy
        print('Calculating occupancy...')
        title = 'Occupancy_'+ time_freq
#        Normalized, S_occupancy =  RFI.spectral_occupancy(freqs,D,outdir,title,3)
        S_occupancy = occupancy(freqs,maskFreq,10**(D/10),2,plot_figs =1)

#        threshold = int(input('Calculate BW loss greater than (percent of the time):  '))
        threshold = (1,2,5,10,30,50,80,90,100)
        for j in threshold:
            BW_loss = np.sum(((S_occupancy>j)))/len(S_occupancy)*100
            print('%0.3f%% of the BW is lost %d%% of the time' % (BW_loss,j))
        
        
    if selection == '7': #reverse Cumulative distrib
        print('Calculating reverse cumulative distrib...')
        P_integ = 10*np.log10(np.sum(10**(D/10),1))
#        P_ave = np.average(P_integ)
#        P_2 = P_integ[P_integ > (P_ave+3)]
        P_2 = P_integ
        fig, ax = plt.subplots(figsize=(8, 4))
        n, bins, patches = ax.hist(P_2, 100, density=True, histtype='step',
                                   cumulative=-1, log=True)
        plt.ylabel('Probability')
        plt.xlabel('Power [dBm]')
        plt.grid()
        title = 'Cumulative distribution function_'+ time_freq
        plt.title(title)
        
        plt.savefig(outdir+ title, dpi=600, bbox_inches='tight')
        
#    os.system('clear')        
    selection = input('\n\nSelect action:\n1: change frequency range\n2: Histogram\n3: Integrated Power\n4: Average\n5: Percentiles\n6: Occupancy\n0: exit\n\nSelect: ')
