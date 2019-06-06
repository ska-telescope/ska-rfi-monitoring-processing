# -*- coding: utf-8 -*-
"""
Created on Wed May 22 02:42:48 2019

@author: F.Divruno
"""
import matplotlib.pyplot as plt
import RFI_general_functions as RFI
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
    
selection = input('\n\nSelect action:\n1: change frequency and time range\n2: Histogram\n3: Integrated Power\n4: Average\n5: Percentiles\n6: Occupancy\n0: exit\n\nSelect: ')
init = 1    
while selection != '0':
    if init: 
        selection = '1'
        init = 0
    if selection == '1': #Change freq range and time
        fmin = float(input('Minimum freq to analyze (MHz): '))
        fmax = float(input('Maximum freq to analyze (MHz): '))
        timestart = input('Start time in sec (enter for 0): ')
        timestop = input('End time in sec (enter for tmax): ')
               
        if timestart == '':
            tmin =  0
        else:
            tmin = int(timestart)
        if timestop == '':
            tmax =  time1[-1]
        else:
            tmax = int(timestop)        

        Daux = data[:,(freq>=fmin) & (freq<=fmax)]/6.5-200 # scaling to match the levels calculated by gianni
        D = Daux[(time1>=tmin) & (time1<=tmax),:]
        Pow = 10*np.log10(np.sum(10**(D/10),1))       
        freqs = freq[(freq>=fmin) & (freq<=fmax)]
        time = time1[(time1>=tmin) & (time1<=tmax)]
        time_freq = str(int((tmax-tmin)/3600))+'hs_'+str(int(fmin)) + 'to'+str(int(fmax))+'MHz'

        print('Done')

    if selection == '2': #Histogram
        print('Calculating histogram of the integrated power in %d to %d MHz' %(fmin,fmax))
        plt.figure()
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
            
        plt.figure()
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
        plt.figure()
        plt.plot(freqs,ave, linewidth=2)
        plt.ylabel('Ampl [dBm]')
        plt.xlabel('freq [MHz]')
        plt.grid()
        title = 'Average_'+ time_freq
        plt.title(title)
        
        plt.savefig(outdir+ title, dpi=600, bbox_inches='tight')


    if selection == '5': #Percentiles
        perc = input('Percentile = ')

        print('Calculating %s percentile' % (perc))
        
        title = 'MRO data '+ time_freq
        perc = int(perc)
        RFI.plot_percentile(freqs,D,perc,outdir,'dBm',title)

 
    if selection == '6': #Occupancy
        print('Calculating occupancy...')
        title = 'Occupancy_'+ time_freq
        S_occupancy =  RFI.spectral_occupancy(freqs,D,outdir,title,1.5)
        threshold = int(input('Calculate BW loss greater than (percent of the time):  '))
        BW_loss = np.sum(((S_occupancy>threshold)))/len(S_occupancy)*100
        print('Loss of %d %% of the BW the %d %% of the time' % (BW_loss,threshold))
        
        
    os.system('clear')        
    selection = input('\n\nSelect action:\n1: change frequency range\n2: Histogram\n3: Integrated Power\n4: Average\n5: Percentiles\n6: Occupancy\n0: exit\n\nSelect: ')
