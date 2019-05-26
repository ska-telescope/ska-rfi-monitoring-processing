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
    
    data = np.zeros([0,29801]).astype('float32')
    N_files = np.size(files)
   
 
    i = 0
    for f in files:
#    for i in range(4): #for debugging, only one file present
#        f = files[i] #for debugging
        fullpath = os.path.join(folder, f)
        if os.path.splitext(fullpath)[1] == '.gz':
           with fits.open(fullpath) as hdul:
               i+=1
               print(str(i)+' of '+str(N_files) + ' Fits files')
               #               hdul.info()
               N = np.size(hdul)
#               for k in range(N):
               for k in range(2): #for debugging
                   try:
                       print(str(k)+' of '+str(N) + ' lines')
                       aux = hdul[k].data
                       data = np.concatenate((data,np.reshape(aux['Amplitude'],[1,29801])),0) #gets the data matrix.
                       
                   except:
                       A=1
    freq = aux['Frequency']    
    return [freq,data]

