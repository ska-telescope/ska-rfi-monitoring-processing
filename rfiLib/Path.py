# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 00:01:50 2019

@author: f.divruno
"""
import numpy as np
import astropy.constants as const
import pyproj as pyproj
import astropy.units as u

ms = 1e-3
us = 1e-6
MHz = 1e6
GHz = 1e9
km = 1e3
minute = 60
hr = 60*minute
km_h = km/hr
k_bolt = 1.38e-23


ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

'''----------------------------------

    ----------------------------------
'''


def Path(posRx,posTx,fc_Hz):
    # as an approximation considers that the emitter is stationary in time
    # in the future movement of the source should be included.
    posRx = np.array(posRx)
    posTx = np.array(posTx)
    
    #range:
    R = np.linalg.norm(posRx - posTx) # R is in km
    R_m = R*km
    FSPL = 20*np.log10(R_m) + 20*np.log10(fc_Hz/MHz) - 27.55 
    
    
    #delay:
    #the delay changes with the angle of incidence and the height of the transmiter.
    lonTx,latTx,altTx = pyproj.transform(ecef,lla,posTx[0]*1e3,posTx[1]*1e3,posTx[2]*1e3) #plot all sats    
    lonRx,latRx,altRx = pyproj.transform(ecef,lla,posRx[0]*1e3,posRx[1]*1e3,posRx[2]*1e3) #plot all sats    
    
    
    #calculate the point on the earth where the phase reference is zero:
    posTxPhaseRef = np.array(pyproj.transform(lla,ecef,lonTx,latTx,altRx))/km
    distRxPhaseRef = np.linalg.norm(posTxPhaseRef-posRx)
    
    #calculate the angle of incidence
    D = posTx-posRx
    alpha0 = (np.arccos(np.dot(posRx,D)/np.linalg.norm(D)/np.linalg.norm(posRx))*180/np.pi)*u.deg
    alpha = 90*u.deg-alpha0
    
    #calculate the distance with respect to the phase centre
    d = distRxPhaseRef*np.cos(alpha)*u.km
    time_delay = (d/const.c).to(u.s).value
    
    return FSPL,time_delay # FSPL in dB , delay in seconds
        
