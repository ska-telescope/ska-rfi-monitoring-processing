# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 00:04:14 2019
    Collection of functions 

@author: f.divruno
"""

import scipy.io as sio
import astropy.coordinates as Coord
import numpy as np

def Coord_to_nparray(coord):
    
    A = np.array([coord[0].si.value,coord[1].si.value,coord[2].si.value])
    
    return A

def save_data(Telescope_list, filename):
    """ Saves data to MATLAB file with vectors named 'time' & 'signal' """
    for i in range(len(Telescope_list)):
#        sio.savemat(filename+'_'+str(i), {"time":Telescope_list[i].time,"ADC_output_rx":Telescope_list[i].ADC_output_rx,"ADC_output_sky":Telescope_list[i].ADC_output_sky })
        # instead of saving the time vector save the sample rate, the time vector can be calculated.
        sio.savemat(filename+'_'+str(i), {"SampleRate":Telescope_list[i].SampleRate,"ADC_output_rx":Telescope_list[i].ADC_output_rx,"ADC_output_sky":Telescope_list[i].ADC_output_sky })
        
def load_data(Telescope_list,*filenames):
    """ Loads data from one or more files, following a generator pattern.
        @yield: time,signal
    """
    for filename in filenames:
        data = sio.loadmat(filename, squeeze_me=True)
        yield (data["time"],data["signal"])
