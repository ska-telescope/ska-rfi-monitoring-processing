# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 23:49:31 2019

@author: f.divruno
"""
import numpy as np
import astropy.coordinates as Coord
import scipy.constants as const

import General as General

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
      
class Sky():
    """
    A class to define sky sources.
    """
    def __init__(self,Name,Src_position,SamplingRate,Duration,Temperature,random_seed = []):
        """            
            Src_position: dict(lon=...*u.deg,lat=...*u.deg) to prove the functionality lon.lat
                          coordinates are given from the source.
            SamplingRate:
            Duration:
            Temperature:
            random_seed:
                
        """
        self.Duration = Duration
        self.SamplingRate = SamplingRate
        self.Temp = Temperature
        self.random_seed = random_seed
        self.BW = SamplingRate
        self.Power = k_bolt*Temperature*self.BW
        self.Ampl = np.sqrt(self.Power/50)
        self.Name = Name
               
        self.lon = Src_position['lon']
        self.lon = Src_position['lat']
        self.Pos = Coord.EarthLocation.from_geodetic(Src_position['lon'],Src_position['lat'])
        
        self.time, self.data = self.sky_signal()
        self.data = self.data*self.Ampl
 
        
        
    def sky_signal(self):
        # Original White noise signal to be summed to every antenna    
        # Sky is generated with an extra time to be able to delay the signal
        # at a maximum angle of 30 deg for a baseline of 80km
        max_baseline = 80*km
        max_angle = 30*np.pi/180
        max_delay = max_baseline/const.c*np.sin(max_angle)
        t_max = self.Duration+2*max_delay
        
        #number of samples
        N = int(t_max*self.SamplingRate)
        time = np.linspace(0,self.Duration,N)
        #generate the random signal
        try:
            np.random.seed(self.random_seed)
            sky = np.random.randn(N)*self.Ampl
        except:
            np.random.seed()
            sky = np.random.randn(N)*self.Ampl
        
        return time,sky
        
    
    # Delay of the star signal depending on the baseline and the angle of the source    
    def rx_sky(self,SamplingRate,duration,pos_rx):
        c = const.c
        pos_src = General.Coord_to_nparray(self.Pos.geocentric)
        cos_phi = np.dot(pos_src,pos_rx)/np.linalg.norm(pos_src)/np.linalg.norm(pos_rx)
        L = np.linalg.norm(pos_src - pos_rx) # R is in m
        L_m = L
        
        
        angle = np.arccos(cos_phi)
        dt = L_m*np.sin(angle)/c
        
        N_delay = int(dt*SamplingRate)
        N = int(duration*SamplingRate)
        
        sky_delayed = self.data[N_delay:(N_delay+N)]
        
        return sky_delayed
