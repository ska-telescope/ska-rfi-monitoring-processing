# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 14:07:34 2019

@author: f.divruno
"""
import os, os.path
import numpy as np
from astropy.io import fits

def read_SKALA2_Bighorns_data(folder):
    files = os.listdir(folder)
    
    data = np.zeros([0,4096]).astype('float32')
    N_files = np.size(files)
   
    i = 0
    for f in files:
#            for f in range(1): #for debugging
#                f = files[0] #for debugging
        fullpath = os.path.join(folder, f)
        if os.path.splitext(fullpath)[1] == '.gz':
           with fits.open(fullpath) as hdul:
               #hdul.info()
               try:
                   data = np.concatenate((data,hdul[0].data),0) #gets the data matrix.
                   i+=1
                   print(str(i)+' of '+str(N_files))
               except:
                   A=1
    
    fullpath = os.path.join(folder, files[0])
    hdul = fits.open(fullpath)  
    N1 = int(hdul[0].header['NAXIS1']) #frequency axis number of points
    N2 = int(hdul[0].header['NAXIS2'])  #Time axis number of points
    fstart = float(hdul[0].header['CRPIX1']) 
    fstep = float(hdul[0].header['CDELT1'])
    freq = np.linspace(fstart,fstart+fstep*N1,N1)
    return [freq,data]
