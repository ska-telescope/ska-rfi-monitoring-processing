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
from TIQ_process_functions import *
from astropy.io import fits


font = {'family' : 'DejaVu Sans','weight' : 'normal','size'   : 22}
matplotlib.rc('font', **font)


#%% Plot the data as waterfall
def plot_waterfall(time,freq,data,title):
      
    fig = plt.figure(figsize=(20,12))
    ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
    #cbax = fig.add_axes((0.85, 0.05, 0.95, .95))
    
    
    left = freq[0]
    right = freq[-1]
    bottom = time[0]
    top = time[-1]
    data_log = (10*np.log10(data))
    cim = ax.imshow( data_log,
        origin='lower', interpolation='nearest',
        cmap= 'jet',
        extent=(left, right, bottom, top),
        )
    
    
    ax.set_aspect(abs(right-left) / abs(top-bottom))
    plt.xlabel('MHz')
    plt.ylabel('time')
    plt.title(title)     


def read_bighorns_SKALA2_data(root):
    files = os.listdir(root)
    
    data = np.zeros([0,4096]).astype('float32')
    N_files = np.size(files)
   
    i = 0
    for f in files:
#            for f in range(1): #for debugging
#                f = files[0] #for debugging
        fullpath = os.path.join(root, f)
        if os.path.splitext(fullpath)[1] == '.gz':
           with fits.open(fullpath) as hdul:
               #hdul.info()
               try:
                   data = np.concatenate((data,hdul[0].data),0) #gets the data matrix.
                   i+=1
                   print(str(i)+' of '+str(N_files))
               except:
                   A=1
    
    fullpath = os.path.join(root, files[0])
    hdul = fits.open(fullpath)  
    N1 = int(hdul[0].header['NAXIS1']) #frequency axis number of points
    N2 = int(hdul[0].header['NAXIS2'])  #Time axis number of points
    fstart = float(hdul[0].header['CRPIX1']) 
    fstep = float(hdul[0].header['CDELT1'])
    freq = np.linspace(fstart,fstart+fstep*N1,N1)
    return [freq,data]

def total_power(time,freq,data,fo,B,Tlow,Thigh,plot_flag):
#    fo = 160 #MHz
#    B = 1 #MHz
#    low_th = 0.4e1
#    high_th = 2e17

    fmin = fo-B/2
    fmax = fo+B/2
    
    fstep = freq[1] - freq[0]
    
    ind1 = int((fmin-freq[0])/fstep)
    ind2 = int((fmax-freq[0])/fstep)
    
    total_power = np.sum(data[:,ind1:ind2]**2,1)
    
    mask = np.ones(len(total_power), dtype=bool)
    for i in range(np.size(total_power)):
        if (total_power[i] < Tlow or total_power[i] > Thigh):
            mask[i] = False
          
    Data2 = total_power[mask]
    time = time[mask]
    
    if plot_flag ==1:
        plt.figure()
        plt.plot(10*np.log10((Data2)))
        plt.title('freq = ' + str(fo) + ' MHz, B = ' + str(B) + ' MHz')
        plt.grid(True,'both')
    return time,Data2
             
    
#%% 


directory = 'C:\\RFI Raw Data\\RFI measurements in AUS\\SKALA 2 with BIGHORNS\\SKALA_Bighorns\\skala_x\\'
 
#get the directories:
dirs = os.listdir(directory)
datenum = 12
[freq,SKALA2_pow] = read_bighorns_SKALA2_data(directory+dirs[datenum]) # the day is the number.

SKALA2_freq = freq
#%% maximum values

plt.figure()
plt.plot(freq,10*np.log10(np.max(SKALA2_pow,0)))
plt.grid(True,'both')
plt.title('SKALA2 with bighorns digitizer')

#%% percentiles
perc = 99
plt.figure()
plt.plot(SKALA2_freq,10*np.log10(np.percentile(SKALA2_pow,perc,axis=0)))
plt.grid(True,'both')
plt.title(str(perc) + ' percentile - SKALA2 with bighorns digitizer')




#%% Waterfall plot
SKALA2_time = np.linspace(0,np.size(SKALA2_pow,0),np.size(SKALA2_pow,0)) # UTC time
SKALA2_time = SKALA2_time + 8*60*60
SKALA2_time[SKALA2_time>24*60*60] = SKALA2_time[SKALA2_time>24*60*60] - 24*60*60
#plot_waterfall(SKALA2_time/3600,freq,SKALA2_pow,'Spectrogram - SKALA2 with bighorns')

#%% total power in a band
# This allows to see what is the time domain behaviour of the interference in this band.

fo = 200 #MHz  # Aeronautic Comms band
B = 300 #MHz

#fo = 128 #MHz  # Aeronautic Comms band
#B = 16 #MHz

#fo = 138.5 #MHz #This is Orbcom
#B = 2 #MHz


#fo = 148.5 #MHz VHF radio?
#B = 1 #MHz

#fo = 255 #MHz Milstar satellite
#B = 30 #MHz

fo = 112.5 #MHz Milstar satellite
B = 9 #MHz


low_th = 0
high_th = 1e50

time2,Data2 = total_power(SKALA2_time*3600,SKALA2_freq,SKALA2_pow,fo,B,low_th,high_th,0)



plt.figure()
plt.plot(time2/3600,10*np.log10(Data2))
plt.grid(True,'both') 
plt.title('SKALA2 with bighorns - ' + dirs[datenum] + ' - Fo=' + str(fo) + ' MHz - B=' + str(B) + 'MHz')
#plt.xlim([0,24])
#plt.ylim([140,220])

#plt.figure()
#plt.hist(10*np.log10(Data2),100)
#plt.grid(True,'both')




#%%
#pop_ind = []
## remove any non wanted date
#for s in dirs:
#    if s[0:4] != '2016':
#        dirs.remove(s)
#    
#    if s[0:8] == '20160908': #remove day 08-09-2016
#        dirs.remove(s)
##        
##    if s[0:8] == '20160925': #remove day 25-09-2016
##        dirs.remove(s)
#        
#datos = np.zeros([0,4096]).astype('float32')
#datos_max = np.zeros([0,4096]).astype('float32')
#datos_98 = np.zeros([0,4096]).astype('float32')
#datos_90 = np.zeros([0,4096]).astype('float32')
#datos_50 = np.zeros([0,4096]).astype('float32')
#datos_min = np.zeros([0,4096]).astype('float32')
#
#for date in dirs:
##for a in range(1): #for debug
##    date = dirs[0] #for debug   
##    try:
#    load_saved = 1
#    if load_saved == 1:
#        try:
#         
#            Dmin = np.load(directory+ '\\'+'SKALA2_1hr_'+date+'.npz')['datos_min']
#            Dmax = np.load(directory+ '\\'+'SKALA2_1hr_'+date+'.npz')['datos_max']
#            D98 = np.load(directory+ '\\'+'SKALA2_1hr_'+date+'.npz')['datos_98']
#            D90 = np.load(directory+ '\\'+'SKALA2_1hr_'+date+'.npz')['datos_90']
#            D50 = np.load(directory+ '\\'+'SKALA2_1hr_'+date+'.npz')['datos_50']
#            #datos = np.load(directory+ '\\'+'SKALA2_2016_'+date+'.npz')['datos']
#            freq = np.load(directory+ '\\'+'SKALA2_1hr_'+date+'.npz')['freq']
#    
#            plt.figure(figsize=(20,12))
#    
#            plt.subplot(221)
#            plt.plot(freq,np.max(D98,0),label='98')
#            plt.ylim([30,150])
#            plt.grid('on')
#            plt.title('98 percentile')
#    
#            plt.subplot(222)
#            plt.plot(freq,np.max(D90,0),label='90')
#            plt.ylim([30,150])
#            plt.grid('on')
#            plt.title('90 percentile')
#            
#            plt.subplot(223)
#            plt.plot(freq,np.max(D50,0),label='50')
#            plt.ylim([30,150])
#            plt.grid('on')
#            plt.title('50 percentile')
#            
#            plt.subplot(224)
#            plt.plot(freq,np.max(Dmin,0),label='max')
#            plt.ylim([30,150])
#            plt.grid('on')
#            plt.title('Min' )
#            
#            plt.suptitle('SKALA2 measurements in 2016 day: ' + date + ' - percentiles over 16min obs')
#    
#            
#            datos_max = np.concatenate((datos_max,Dmax),0)
#            datos_min = np.concatenate((datos_min,Dmin),0)
#            datos_98 = np.concatenate((datos_98,D98),0)
#            datos_90 = np.concatenate((datos_90,D90),0)
#            datos_50 = np.concatenate((datos_50,D50),0)
#        except:
#            print('File ' + date + ' not found')
#    else:
#    
#        root = directory + '\\' + date
#        try: 
#            files = os.listdir(root)
#    # Each file is a matrix of 1000 seconds, 4096 freq points, thats 16 minutes aprox
#    # there are aprox 4 files per hour.
#    # Data reduction:
#    # In each file look for max, min, 98%, 90%, 50%, limt line for RFI.
#    # Save the vectors for these percetiles and then plot the maximum levels for all the data capture available.
#    #
#            fstart = 2
#            #freq = np.linspace(fstart,fstart+fstep*N1,int(N1))
#        
#            datos_t = np.zeros([0,4096]).astype('float32')
#            datos_max_t = np.zeros([0,4096]).astype('float32')
#            datos_98_t = np.zeros([0,4096]).astype('float32')
#            datos_90_t = np.zeros([0,4096]).astype('float32')
#            datos_50_t = np.zeros([0,4096]).astype('float32')
#            datos_min_t = np.zeros([0,4096]).astype('float32')
#            aux_t = np.zeros([0,4096]).astype('float32')
#            
#            
#        #    plt.figure()
#            N_frames = 4
#            i = 0
#            for f in files:
##            for f in range(1): #for debugging
##                f = files[0]
#                if i == N_frames: # toma 4 archivos de 16 minutos  para hacer una captura de tiempo equivalente a 1 hora.
#                    aux_t = np.zeros([0,4096]).astype('float32')
#                    i = 0
#                    
#                fullpath = os.path.join(root, f)
#                if os.path.splitext(fullpath)[1] == '.gz':
#                   with fits.open(fullpath) as hdul:
#                       hdul.info()
#                       try:
#                           aux_t = np.concatenate((aux_t,hdul[0].data),0) #gets the data matrix.
#                           i+=1
#                           
#                                                      
#                           if i == N_frames:
#                               aux_sort = 10*np.log10(np.sort(aux_t,0)) # pasa a logaritmico
#                
#                               index = int(np.size(aux_sort,0)*0/100)
#                               aux_min = aux_sort[index,:]
#                
#                               index = np.size(aux_sort,0)-1
#                               aux_max = aux_sort[index,:]
#                
#                               index = int(np.size(aux_sort,0)*98/100)
#                               aux_98 = aux_sort[index,:]
#                
#                               index = int(np.size(aux_sort,0)*90/100)
#                               aux_90 = aux_sort[index,:]
#                
#                               index = int(np.size(aux_sort,0)*50/100)
#                               aux_50 = aux_sort[index,:]
#                               
#                                       
#                               datos_max_t = np.concatenate((datos_max_t,np.reshape(aux_max,[1,4096])),0)
#                               datos_min_t = np.concatenate((datos_min_t,np.reshape(aux_min,[1,4096])),0)
#                               datos_98_t = np.concatenate((datos_98_t,np.reshape(aux_98,[1,4096])),0)
#                               datos_90_t = np.concatenate((datos_90_t,np.reshape(aux_90,[1,4096])),0)
#                               datos_50_t = np.concatenate((datos_50_t,np.reshape(aux_50,[1,4096])),0)
#                               
#    #                           datos = np.concatenate((datos,10*np.log10(aux)),0)
#                                   
#                               N1 = int(hdul[0].header['NAXIS1']) #frequency axis number of points
#                               N2 = int(hdul[0].header['NAXIS2'])  #Time axis number of points
#                               fstart = float(hdul[0].header['CRPIX1']) 
#                               fstep = float(hdul[0].header['CDELT1'])
#                               freq = np.linspace(fstart,fstart+fstep*N1,N1)
#                               
##                               plot_waterfall(freq,aux_t)
#                               
#                       except:
#                           A=1    
#        
#                       
#            np.savez_compressed(directory + '\\'+'SKALA2_1hr_' + date, freq=freq, datos_max=datos_max_t, datos_min=datos_min_t,datos_50=datos_50_t, datos_90=datos_90_t,datos_98=datos_98_t)
#    
#        except:
#            A=1
#            
##plt.figure()
##plt.plot(freq,np.transpose(datos))
##
##plt.figure()
##plt.plot(freq,np.transpose(aux_sort))
##    
##plt.figure()
##plt.plot(freq,np.transpose(baseline))
#
#
##%% 
#plt.figure(figsize=(20,12))
#plt.plot(freq,np.max(datos_max,0),label='max')
#plt.plot(freq,np.max(datos_98,0),label='98')
#plt.plot(freq,np.max(datos_90,0),label='90')
#plt.plot(freq,np.max(datos_50,0),label='50')
#plt.plot(freq,np.max(datos_min,0),label='min')
#plt.ylim([30,150])
#plt.legend()
#plt.grid('on')
#plt.title('')
#
#
#plt.figure(figsize=(20,12))
#plt.plot(freq,np.max(datos_98,0),label='98')
#plt.ylim([30,150])
#plt.legend()
#plt.grid('on')
#plt.title('SKALA2 - 2016 : Maximum levels of percentile in 16min captures ' + date)
#
#plt.figure(figsize=(20,12))
#plt.plot(freq,np.max(datos_90,0),label='90')
#plt.ylim([30,150])
#plt.legend()
#plt.grid('on')
#plt.title('SKALA2 - 2016 : Maximum levels of percentile in 16min captures' + date)
#
#
#plt.figure(figsize=(20,12))
#plt.plot(freq,np.max(datos_50,0),label='50')
#plt.ylim([30,150])
#plt.legend()
#plt.grid('on')
#plt.title('SKALA2 - 2016 : Maximum levels of percentile in 16min captures' + date)
#
#
#
##plt.figure()
##plt.plot(freq,aux_max_d,label='max')
##plt.plot(freq,aux_98_d,label='98')
##plt.plot(freq,aux_90_d,label='90')
##plt.plot(freq,aux_50_d,label='50')
##plt.legend()
##plt.title('differences with base line')
#
#
#
          