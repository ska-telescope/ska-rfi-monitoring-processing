# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 00:01:50 2019

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
'''----------------------------------

    ----------------------------------
'''


def Path(pos_rx,pos_tx,fc_Hz):
    # as an approximation considers that the emitter is stationary in time
    # in the future movement of the source should be included.
    
    v_light = 3e8 #approximation of light speed.
    
    R = np.linalg.norm(pos_rx - pos_tx) # R is in km
    R_m = R*km
    FSPL = 20*np.log10(R_m) + 20*np.log10(fc_Hz/MHz) - 27.55 
    time_delay = R_m/v_light
    
    return FSPL,time_delay # FSPL in dB , delay in seconds
        
