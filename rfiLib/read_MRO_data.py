# -*- coding: utf-8 -*-
"""
Created on Mon May  6 06:19:33 2019

@author: f.divruno
"""

import os, os.path
import numpy as np
from astropy.io import fits

def read_MRO_data(indir,outdir,ext='.fits'):
    # use ext = '.fits' if the files are decompressed
    # use ext = '.gz' if files are compressed (significantly more time)
    # this function does read ALL the files inside the input directory.
    
    files = os.listdir(indir)
    freq_points = 29801
    data = np.zeros([0,freq_points]).astype('float32')
#    data_file = np.zeros([0,freq_points]).astype('float32')
    data_file = np.zeros([680,freq_points]).astype('float32')
    N_files = np.size(files)
   
 
    i = 0
    for f in files:
#    for i in range(2): #for debugging
#        f = files[i] #for debugging
        fullpath = os.path.join(indir, f)
        if os.path.splitext(fullpath)[1] == ext:            
           with fits.open(fullpath) as hdul:
               i+=1
               print(str(i)+' of '+str(N_files) + ' Fits files')
               #               hdul.info()
               N = np.size(hdul)
               j = 0
               for k in range(N):
#               for k in range(2): #for debugging
                   try:
                       print(str(k)+' of '+str(N) + ' lines')
#                       aux = hdul[k].data
                       data_file[j,:] = hdul[k].data['Amplitude']
                       j+=1
#                       data_file = np.concatenate((data_file,np.reshape(aux['Amplitude'],[1,freq_points])),0) #gets the data matrix.
                       
                   except:
                       A=1
               freq = hdul[k].data['Frequency']
               print('Saving the file as npz...')
#               np.savez_compressed(outdir + 'MRO_rfidata_' + str(i), freq=freq, data_file=data_file)
               data = np.concatenate((data,data_file),0)
#               data_file = np.zeros([0,freq_points]).astype('float32')
               
                       
    return [freq,data]

