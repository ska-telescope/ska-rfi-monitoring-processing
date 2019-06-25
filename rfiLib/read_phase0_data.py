# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 12:20:01 2019

@author: f.divruno
"""
import numpy as np
import os, os.path

def read_phase0_data(directory,files_num):
    # files_num provide the number of files to read inside the supplied directory.

    files = os.listdir(directory)
        
    header = 131072
    if files_num == 'all':
        N_files = np.size(files)
    else:
        N_files = files_num
    time_vect = np.zeros([N_files])
    data = np.zeros([N_files,header]).astype('int8')
    i = 0
    for i in range(N_files):
        f = files[i]
        fullpath = os.path.join(directory, f)
        if os.path.splitext(fullpath)[1] == '.tdd':
           with open(fullpath,'r') as f:
               header = np.fromfile(f, dtype=np.dtype('>f8') , count = 1)#dtype=np.floating point 8 bytes
               data[i,:] = np.fromfile(f, dtype=np.dtype('>i1') , count = 131072)#dtype=np.int64
        sec = os.path.splitext(fullpath)[0][-2::]
        mins = os.path.splitext(fullpath)[0][-4:-2]
        hours = os.path.splitext(fullpath)[0][-6:-4]
        time_vect[i] = int(hours)*60*60+int(mins)*60+int(sec)
        i += 1
        print(str(i))
        
    return time_vect,data
