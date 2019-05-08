# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 12:21:22 2018
Open and process FITS files produced by BIGHORNS system

@author: f.divruno
"""


import os, os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#from TIQ_process_functions import *
import RFI_general_functions as RFI
#from astropy.io import fits


font = {'family' : 'DejaVu Sans','weight' : 'normal','size'   : 22}
matplotlib.rc('font', **font)


#%% SKALA2 with bighorns
#
#
#directory = 'C:\\RFI Raw Data\\RFI measurements in AUS\\SKALA 2 with BIGHORNS\\SKALA_Bighorns\\skala_x\\'
# 
##get the directories:
#dirs = os.listdir(directory)
#datenum = 1
#[SKALA2_freq,SKALA2_pow] = RFI.read_SKALA2_Bighorns_data(directory+dirs[datenum]) # the day is the number.

#%% MRO Data ---- testing
filepath  = 'C:\\RFI Raw Data\\RFI measurements in AUS\\MRO_monitoring\\'
[MOR_freq,MRO_E] = RFI.read_MRO_data(filepath)


#%% maximum values

f= SKALA2_freq
data= SKALA2_pow
title= 'SKALA2 with bighorns digitizer - '+dirs[datenum]
RFI.plot_percentile(f,data,100,title='Maximum values - '+title)



#%% percentiles
perc = 98
f= SKALA2_freq
data= SKALA2_pow
title= 'SKALA2 with bighorns digitizer- '+dirs[datenum]
RFI.plot_percentile(f,data,perc,title)



#%% spectrogram
SKALA2_time = np.linspace(0,np.size(SKALA2_pow,0),np.size(SKALA2_pow,0)) # UTC time
SKALA2_time = SKALA2_time + 8*60*60
SKALA2_time[SKALA2_time>24*60*60] = SKALA2_time[SKALA2_time>24*60*60] - 24*60*60

Fstart = 0
Fstop = 30
time = SKALA2_time
f= SKALA2_freq
data= SKALA2_pow
title= 'SKALA2 with bighorns digitizer- '+dirs[datenum]

RFI.plot_spectrogram(time/3600,f,data,'Spectrogram -'+title,Fstart,Fstop,Tstart=0,Tstop=0)


#%% total power in a band
# This allows to see what is the time domain behaviour of the interference in this band.
time = SKALA2_time
f= SKALA2_freq
data= SKALA2_pow
title= 'SKALA2 with bighorns digitizer'

#fo = 200 #MHz  # Full Band
#B = 300 #MHz

fo = 128 #MHz  # Aeronautic Comms band
B = 16 #MHz

#fo = 138.5 #MHz #This is Orbcom
#B = 2 #MHz


#fo = 148.5 #MHz VHF radio?
#B = 1 #MHz

#fo = 255 #MHz Milstar satellite
#B = 30 #MHz

#fo = 112.5 #MHz Milstar satellite
#B = 9 #MHz

#fo = 16 #MHz  # low freq stuf
#B = 25 #MHz


low_th = 1e9
high_th = 1e12

time2,Data2 = RFI.total_power(time*3600,f,data,fo,B,low_th,high_th,0)

plt.figure()
plt.plot(time2/3600,10*np.log10(Data2))
plt.grid(True,'both') 
plt.title(title + ' - Fo=' + str(fo) + ' MHz - B=' + str(B) + 'MHz')


          