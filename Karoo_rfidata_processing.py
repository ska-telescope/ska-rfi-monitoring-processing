# -*- coding: utf-8 -*-
"""
Created on Wed May 22 02:42:48 2019

@author: F.Divruno
"""
import matplotlib.pyplot as plt
import matplotlib
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

location='M48'
location='SKA021'
location='SKA004'

if ubuntu:
# in ubuntu
    matplotlib.use('Agg')
    plt.close('all')
    outdir = r'/mnt/data/SARAO_2017_meas/'+location+'/results/'
    indir = r'/mnt/data/SARAO_2017_meas/'+location+'/'
else:
# in windows
#    indir = 'C:\\Users\\F.Divruno\\Dropbox (SKA)\\14- RFI environment\\02- ZA\\rfidata\\'
#    outdir = 'C:\\Users\\F.Divruno\\Dropbox (SKA)\\14- RFI environment\\02- ZA\\rfidata\\results\\'
    indir = 'C:\\Users\\F.Divruno\\Dropbox (SKA)\\14- RFI environment\\02- ZA\\rfidata\\'
    outdir = 'C:\\Users\\F.Divruno\\Dropbox (SKA)\\14- RFI environment\\02- ZA\\rfidata\\results\\'

#%% Read the files from disk
    
read_files = int(input('---Input data--\n\nData in memory (0)\nRead saved npz file (1)\nRead HDF52 files (2)\n Selection: '))
if read_files==2:
    print('reading files...')
    [freq,data] = RFI.read_hdf5_data(indir,outdir)
    
    if int(input('Save npz file?\nYes(1) No(0): ')):
        np.savez_compressed(outdir + r'ZA_rfidata_full', freq=freq, data=data)
    
elif read_files==1:
    try:
        print('Loading the complete set of saved data...')
        Aux = np.load(outdir+r'ZA_rfidata_full.npz')
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
t_step = 10
tmax = (len(data)-1)*t_step
time1 = np.linspace(0,tmax,len(data))
    
selection = input('\n\nSelect action:\n1: change frequency and time range\n2: Histogram\n3: Integrated Power\n4: Average\n5: Percentiles\n6: Occupancy\n0: exit\n\nSelect: ')
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
        
        Daux = data[:,i_f_start:i_f_stop] 
        D = Daux[i_t_start:i_t_stop,:]
        print('Data sliced')
        freqs = freq[i_f_start:i_f_stop]
        time = time1[i_t_start:i_t_stop]
        print('freq and time sliced')
        Pow = 10*np.log10(np.sum(10**(D/10),1))       
        time_freq = str(int((tmax-tmin)/3600))+'hs_'+str(int(fmin)) + 'to'+str(int(fmax))+'MHz'

        print('Done')
    if selection == '2': #Histogram
        print('Calculating histogram of the integrated power in %d to %d MHz' %(fmin,fmax))
        plt.figure(figsize=(30,20))
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
            
        plt.figure(figsize=(30,20))
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
        plt.figure(figsize=(30,20))
#        plt.figure()
        plt.plot(freqs,ave, linewidth=2)
        plt.ylabel('Ampl [dBm]')
        plt.xlabel('freq [MHz]')
        plt.grid()
        title = 'Average_'+ time_freq
        plt.title(title)
        
        plt.savefig(outdir+ title, dpi=600)#, bbox_inches='tight')


    if selection == '5': #Percentiles00
        perc = input('Percentile = ')

        print('Calculating %s percentile' % (perc))
        
        title = 'MRO data '+ time_freq
        perc = int(perc)
        RFI.plot_percentile(freqs,D,perc,outdir,'dBm',title)

 
    if selection == '6': #Occupancy
        print('Calculating occupancy...')
        title = 'Occupancy_'+ time_freq
        Normalized, S_occupancy =  RFI.spectral_occupancy(freqs,D,outdir,title,2)
#        threshold = int(input('Calculate BW loss greater than (percent of the time):  '))
        threshold = (1,2,5,10,30,50,80,90,100)
        for j in threshold:
            BW_loss = np.sum(((S_occupancy>j)))/len(S_occupancy)*100
            print('%0.3f%% of the BW is lost %d%% of the time' % (BW_loss,j))
        
        
#    os.system('clear')        
    selection = input('\n\nSelect action:\n1: change frequency range\n2: Histogram\n3: Integrated Power\n4: Average\n5: Percentiles\n6: Occupancy\n0: exit\n\nSelect: ')
