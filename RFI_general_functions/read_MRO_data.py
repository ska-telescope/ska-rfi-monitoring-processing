# -*- coding: utf-8 -*-
"""
Created on Mon May  6 06:19:33 2019

@author: f.divruno
"""

import os, os.path
import numpy as np
from astropy.io import fits

def read_MRO_data(folder):
    files = os.listdir(folder)
    
    data = np.zeros([0,1913601]).astype('float32')
    N_files = np.size(files)
   
    i = 0
    for f in files:
#            for f in range(1): #for debugging
#                f = files[0] #for debugging
        fullpath = os.path.join(folder, f)
        if os.path.splitext(fullpath)[1] == '.gz':
           with fits.open(fullpath) as hdul:
               #hdul.info()
               for k in range(np.size(hdul)):
                   try:
                       aux = hdul[k].data
                       data = np.concatenate((data,np.reshape(aux['Amplitude'],[1,1913601])),0) #gets the data matrix.
                       print(str(k)+' of '+str(np.size(hdul)) + ' lines')
                   except:
                       A=1
               i+=1
               print(str(i)+' of '+str(N_files))
    freq = aux['Frequency']    
    return [freq,data]
