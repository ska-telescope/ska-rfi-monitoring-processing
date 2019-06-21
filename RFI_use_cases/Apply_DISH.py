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
        

def Apply_DISH(Telescope_list,Band='B1',scaling = 'Correlator_opimized', Gain = 0):
    """
        scaling: 
                Correlator_optimized scales the input signals so tht the noise (without RFI) scales to 0.335*Stdv(noise) = 1 level of ADC.
                Linearity_optimized scales the input signals so that the RMS power (noise + RFI) scales to the full scale of the ADC (minimum clipping)
                Defined_Gain scales the input signal by a defined gain, parameter Gain needs to be provided 
        Band: 'B1' 'B2' 'B3' 'B4' 'B5a' 'B5b'
        
    """
    SampleRate = Telescope_list[0].SampleRate
    Duration = Telescope_list[0].Duration
    for i in len(Telescope_list):
        Telescope_list[i].Noise_scale(scaling,Gain) # Calculate the gain that the system needs to have to get an optimum correlator performance.
        
        
    
    