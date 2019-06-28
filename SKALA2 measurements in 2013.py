# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 16:34:12 2018
Processing of SKALA2 measurements in 2013.



@author: f.divruno
"""

import numpy as np
import matplotlib.pyplot as plt
#from TIQ_process_functions import *
import rfiLib as RFI

#%%

root_dir = 'C:\\Users\\f.divruno\\Dropbox (SKA)\\14- RFI environment\\01- Australia\\SKALA2 data 2013\\'
archivs = ['Pos1_Input2_Segment1', ]
freq = np.array([])
power = np.array([])
time = np.array([])
time_aux = np.array([])
power_aux = np.array([])
f1 = np.array([])
f2 = np.array([])
i = 0
percent = np.array([])
freq = np.array([])

margin = 1 #6 dB margin
n2 = 4950


f1 = np.linspace(2,2+2001*498/10000,2001) #this is the first segment of the data measured.
f = np.linspace(2,500,10001) #this is the freq range
freq3 = f

#%% Read the files from the disk
read_files =0
if read_files:
    Pos1_Input1 = np.array(np.zeros([10001,5000])) # saves the space
    Pos1_Input2 = np.array(np.zeros([10001,5000])) # saves the space
    for j in range(5):
        file = 'Pos1_Input1_Segment' + str(j+1)
        with open(root_dir + file) as archivo:
            aux = archivo.readlines()
            k = 0
            print(file)
            for line in aux:
                aux2 = line.split(' ')
                aux2.pop(0)
                aux2[-1] = aux2[-1].split('\n')[0]
                power_aux = np.asarray([float(i) for i in aux2]) #convert a line of numbers in text to numpy array
                Pos1_Input1[k+2000*j,:] = power_aux
                k = k+1
    for j in range(5):
        file = 'Pos1_Input2_Segment' + str(j+1)
        with open(root_dir + file) as archivo:
            aux = archivo.readlines()
            k = 0
            print(file)
            for line in aux:
                aux2 = line.split(' ')
                aux2.pop(0)
                aux2[-1] = aux2[-1].split('\n')[0]
                power_aux = np.asarray([float(i) for i in aux2]) #convert a line of numbers in text to numpy array
                Pos1_Input2[k+2000*j,:] = power_aux
                k = k+1
    
    # save the matrices in disk
    np.savez_compressed(root_dir+'SKALA2_Pos1_Input1', Pos1_Input1 = Pos1_Input1)
    np.savez_compressed(root_dir+'SKALA2_Pos1_Input2', Pos1_Input2 = Pos1_Input2)
else:
    Pos1_Input1 = np.load(root_dir+'SKALA2_Pos1_Input1.npz')['Pos1_Input1']
    Pos1_Input2 = np.load(root_dir+'SKALA2_Pos1_Input2.npz')['Pos1_Input2']


SKALA2_2013_P1_I1 = 10**(np.transpose(Pos1_Input1)/10)
SKALA2_2013_P1_I2 = 10**(np.transpose(Pos1_Input2)/10)

#%% maximum values

RFI.plot_percentile(freq3,SKALA2_2013_P1_I1,100,title='Maximum values - SKALA2 measurements in 2013 - Pos1 Input1')
RFI.plot_percentile(freq3,SKALA2_2013_P1_I2,100,title='Maximum values - SKALA2 measurements in 2013 - Pos1 Input2')


#%% percentiles
RFI.plot_percentile(freq3,SKALA2_2013_P1_I1,98,title=str(perc) + ' percentile - SKALA2 in 2013 - Pos1 Input1')
RFI.plot_percentile(freq3,SKALA2_2013_P1_I1,98,title=str(perc) + ' percentile - SKALA2 in 2013 - Pos1 Input2')


#%% spectrogram

SKALA2_time2 = np.linspace(0,np.size(SKALA2_2013_P1_I1,0)*10*60,np.size(SKALA2_2013_P1_I1,0)) # 10 minutes steps

Fstart = 0
Fstop = 30

RFI.plot_spectrogram(SKALA2_time2/3600,freq3,SKALA2_2013_P1_I1,'Spectrogram - SKALA2 2013',Fstart,Fstop,Tstart=0,Tstop=0)


#%% total power in a band
# This allows to see what is the time domain behaviour of the interference in this band.

fo = 200 #MHz  # full band
B = 398 #MHz

#fo = 128 #MHz  # Aeronautic Comms band
#B = 16 #MHz

#fo = 138.5 #MHz #This is Orbcom
#B = 2 #MHz


#fo = 148.5 #MHz VHF radio?
#B = 1 #MHz

#fo = 255 #MHz Milstar satellite
#B = 30 #MHz

#fo = 112.5 #MHz Milstar satellite
#B = 9 #MHz
#
#fo = 16 #MHz  # low freq stuf
#B = 25 #MHz



low_th = 0
high_th = 1e200

time2,Data2 = RFI.total_power(SKALA2_time2*3600,freq3,SKALA2_2013_P1_I2,fo,B,low_th,high_th,0)

plt.figure()
plt.plot(time2/3600,10*np.log10(Data2))
plt.grid(True,'both') 
plt.title('SKALA2 in 2013 - Pos1 Input1 - Fo=' + str(fo) + ' MHz - B=' + str(B) + 'MHz')

