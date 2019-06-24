# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 15:44:41 2019

@author: f.divruno
"""

import scipy.signal as Sig
import numpy as np
import matplotlib.pyplot as plt
import RFI_general_functions as RFI
import astropy.coordinates as Coord
import astropy.units as u
from poliastro.twobody import Orbit
from poliastro.bodies import Earth
import siggen as siggen
import scipy.constants as const


ms = 1e-3
us = 1e-6
MHz = 1e6
GHz = 1e9
km = 1e3
minute = 60
hr = 60*minute
km_h = km/hr
k_bolt = 1.38e-23
        

def Apply_DISH(Telescope_list,Band='B1',scaling = 'Correlator_opimized', atten = 0):
    """
        scaling: 
                Correlator_optimized scales the input signals so tht the noise (without RFI) scales to 0.335*Stdv(noise) = 1 level of ADC.
                Linearity_optimized scales the input signals so that the RMS power (noise + RFI) scales to the full scale of the ADC (minimum clipping)
                Defined_Gain scales the input signal by a defined atten, parameter atten needs to be provided 
        Band: 'B1' 'B2' 'B3' 'B4' 'B5a' 'B5b'
        
    """
#    SampleRate = Telescope_list[0].SampleRate
#    Duration = Telescope_list[0].Duration

    for i in range(len(Telescope_list)):
        #filter and scales the signals according to the Band and optimiztion selected, attenuation can be provided.
        Telescope_list[i].Apply_analog_chain(Band,scaling,atten=0,f_offset=0) 
        # The signal inputing to the ADC is in the variable Receiver.ADC_input
        
        # digitize the signals.
        Telescope_list[i].Apply_ADC(nBits=12)
        # The output signal is stored in Receiver.ADC_output        
    
    return Telescope_list
    