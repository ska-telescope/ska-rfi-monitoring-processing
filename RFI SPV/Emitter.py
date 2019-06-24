# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 23:49:59 2019

@author: f.divruno
"""
import numpy as np
import astropy.coordinates as Coord

from Signal import Signal

ms = 1e-3
us = 1e-6
MHz = 1e6
GHz = 1e9
km = 1e3
minute = 60
hr = 60*minute
km_h = km/hr
k_bolt = 1.38e-23

class Emitter():
    #Generate emitter object
    
    def __init__(self,Name,Emit_type,Pos_ini, duration,SampleRate,position_params):
        #Initializaition:
        self.Duration = duration
        self.Name = Name
        self.Emit_type = Emit_type
        self.SampleRate = SampleRate
        self.position_params = position_params 
        self.Pos_ini = Pos_ini
        
        if Emit_type == "Aeroplane":
            self.Signals = list()
            DME_power = 60 #dBm = 1kW
            ADSB_power = 53 #dBm = 250W
            self.Signals.append(Signal('DME',self.Duration,self.SampleRate,1100*MHz,DME_power,0,random_seed=[]))
            self.Signals.append(Signal('ADS-B',self.Duration,self.SampleRate,1090*MHz,ADSB_power,0,random_seed=[]))
        
        if Emit_type == 'Sky':
            self.Signals = list()
            scale = 1 # needs some scaling that makes sense 
            self.Signals.append(scale*self.sky_signal(self.SampleRate,self.Duration))
                        
            
        self.Atitude()

    
    def Atitude(self):
        # generates (x,y,z) coordinates in time
        if self.Emit_type == "Aeroplane":
            self.Pos = self.propagate_aeroplane()

        if self.Emit_type == "Sky":
            self.Pos = self.propagate_aeroplane()

    
#        if Emit_type == "NGSO_Sat":
##            [x,y,z] = propagate_NGSO()
##            
#        if Emit_type == "GSO_Sat":
##            [x,y,z] = propagate_GSO()
##            
#        if Emit_type == "MSO_Sat":
##            [x,y,z] = propagate_MSO()
##            
#        if Emit_type == "Ground_mobile":    
##            [x,y,z] = propagate_GNDmobile()
##        
#        if Emit_type == "Ground":    
##            [x,y,z] = init_params['xyz']


    def propagate_aeroplane(self):
        
        self.time = np.linspace(0,self.Duration) #50 samples of position.
        
        height_i = self.Pos_ini['height_i']
        lat_i = self.Pos_ini['lat_i']
        lon_i = self.Pos_ini['lon_i']
        
#        height_i = 10*u.km
#        lat_i = -30*u.deg
#        lon_i = 20*u.deg
        
        Pos_ini = Coord.EarthLocation.from_geodetic(lon_i,lat_i,height_i)
        
        x = np.ones(len(self.time))*Pos_ini.geocentric[0]
        y = np.ones(len(self.time))*Pos_ini.geocentric[1]
        z = np.ones(len(self.time))*Pos_ini.geocentric[2]
        
        return list([x,y,z])
    
  