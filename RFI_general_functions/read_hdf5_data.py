# -*- coding: utf-8 -*-
"""
Created on Mon May  6 06:19:33 2019

@author: f.divruno
"""

import os, os.path
import numpy as np
#from astropy.io import fits
import h5py as h5py


def read_hdf5_data(indir,outdir,ext='.h5'):
    # use ext = '.fits' if the files are decompressed
    # use ext = '.gz' if files are compressed (significantly more time)
    # this function does read ALL the files inside the input directory.
    
    files = os.listdir(indir)
    N_files = np.size(files)
   
 
    i = 0
    for f in files:
#    for i in range(2): #for debugging
#        f = files[i] #for debugging
        fullpath = os.path.join(indir, f)
        if os.path.splitext(fullpath)[1] == ext:            
           with h5py.File(fullpath) as HDF:
               i+=1
               print(str(i)+' of '+str(N_files) +' ' + ext + ' files')
               freq = list(HDF['freqs'])
               data = list(HDF['calibrated_spectrum'])
               
               dat = np.ndarray([len(data),len(data[0])])
               for j in range(len(data)):
                   dat[j] = data[j]
        freqs = np.asarray(freq)/1e6 # to MHz
#        np.savez_compressed(outdir + 'ZA_rfidata_' + str(i), freq=freq, data_file=data_file)
                       
    return [freqs,dat]

