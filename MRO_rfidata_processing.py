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
fmin = float(input('Minimum freq to analyze (MHz): '))
fmax = float(input('Maximum freq to analyze (MHz): '))
D = data[:,(freq>=fmin) & (freq<=fmax)]/6.5-200 # scaling to match the levels calculated by gianni
freqs = freq[(freq>=fmin) & (freq<=fmax)]
Pow = 10*np.log10(np.sum(10**(D/10),1))


#%% Calculate histogram

if input('Calculate histogram? (enter=no)') != '':
    print('Calculating histogram')
    plt.figure()
    plt.hist(Pow,500)
    plt.xlabel('Ampl [dBm]')
    plt.savefig(outdir+ 'power_histogram_freq_'+str(int(fmin)) + 'to'+str(int(fmax))+'' , dpi=500, bbox_inches='tight')

if input('Calculate integrated power? (enter=no)') != '':
    timestart = input('Start time (enter for 0): ')
    timestop = input('Start time (enter for tmax, x to continue): ')
    t_step = 1200/680
    tmax = (len(Pow)-1)*t_step
    time = np.linspace(0,tmax,len(Pow))
        
    while  timestop != 'x':
        print('Calculating total power...')
        
        if timestart == '':
            tmin =  0
        else:
            tmin = int(timestart)

        if timestop == '':
            tmax =  time[-1]
        else:
            tmax = int(timestop)
        
        time2 = time[(time>tmin) & (time<=tmax)]
        Pow2 = Pow[(time>tmin) & (time<=tmax)]
            
        plt.figure()
        plt.plot(time2,Pow2)
        plt.xlabel('time [seconds]')
        plt.ylabel('Ampl [dBm]')
        title = 'total_power_'+str(tmin)+'sec_to_'+str(int(tmax))+'sec_'+str(int(fmin)) + 'to'+str(int(fmax))+'MHz'
        plt.title(title)
        plt.savefig(outdir+title , dpi=500, bbox_inches='tight')
        print('Done')
        timestart = input('Start time (enter for 0): ')
        timestop = input('Start time (enter for tmax, x to continue): ')


if input('Calculate Average? (enter=no)') != '':
    print('Calculating Average')
    ave = np.average(D,0)            
    plt.figure()
    plt.plot(freqs,ave)
    plt.ylabel('Ampl [dBm]')
    plt.xlabel('freq [MHz]')
    plt.savefig(outdir+ 'Average_all', dpi=100, bbox_inches='tight')


if input('Calculate percentile? (enter=no)') != '':
    perc = input('Percentile? (x to exit): ')
    while  perc != 'x':
        print('Calculating %s percentile' % (perc))
        
        title = 'MRO data '+str(int(fmin)) + ' to '+str(int(fmax)) + ' MHz'
        perc = int(perc)
        RFI.plot_percentile(freqs,D,perc,outdir,'dBm',title)
        perc = input('Percentile? (x to exit): ')
    
#
#title = 'MRO data'
#perc = 98
#RFI.plot_percentile(freqs,D,perc,outdir,'dBm',title)
#
#
#title = 'MRO data'
#perc = 95
#RFI.plot_percentile(freqs,D,perc,outdir,'dBm',title)
#
#
#title = 'MRO data'
#perc = 90
#RFI.plot_percentile(freqs,D,perc,outdir,'dBm',title)
#
#title = 'MRO data'
#perc = 50
#RFI.plot_percentile(freqs,D,perc,outdir,'dBm',title)
#

#%% Calculate spectral occupancy

#S_occupancy =  RFI.spectral_occupancy(freqs,D,outdir,1.5)

