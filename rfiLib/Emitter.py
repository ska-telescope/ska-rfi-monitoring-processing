# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 23:49:59 2019

@author: f.divruno
"""
import numpy as np
import astropy.coordinates as Coord

from rfiLib.Signal import Signal

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
    
    def __init__(self,Name,Emit_type,Pos_ini, duration,SampleRate,position_params,random_seed=[]):
        #Initializaition:
        self.Duration = duration
        self.Name = Name
        self.Emit_type = Emit_type
        self.SampleRate = SampleRate
        self.position_params = position_params 
        self.Pos_ini = Pos_ini
        self.randomSeed = random_seed
        
        if Emit_type == "Airplane":
            self.Signals = list()
            DME_power = 60 #dBm = 1kW
            ADSB_power = 53 #dBm = 250W
            np.random.seed(int(self.randomSeed*100)) #loads the seed in case of wanting to repeat the same signal.
            freqDme = (np.random.random(1)-0.5)*62.5 + 1087.5 #DME frequencies for tx airplane : 1025-1150MHz
            self.Signals.append(Signal('DME',self.Duration,self.SampleRate,freqDme*MHz,DME_power,0,self.randomSeed))
            self.Signals.append(Signal('ADS-B',self.Duration,self.SampleRate,1090*MHz,ADSB_power,0,self.randomSeed*5))
        elif Emit_type == 'Sky':
            self.Signals = list()
            scale = 1 # TO-DO :  needs some scaling that makes sense 
            self.Signals.append(scale*self.sky_signal(self.SampleRate,self.Duration))
        else:
            raise Exception("Can't instantiate emitter. Illegal RFI emitter type: " + Emit_type)                
            
        self.Atitude()

    
    def Atitude(self):
        # generates (x,y,z) coordinates in time
        if self.Emit_type == "Airplane":
            self.Pos = self.propagate_aeroplane()
        elif self.Emit_type == "Sky":
            self.Pos = self.propagate_aeroplane()
        else:
            raise Exception("Can't set Atitude position. Illegal RFI emitter type: " + self.Emit_type)                

    
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
        
        # TO-DO make some propagation function that can move the emitter, not done yet.
        
        x = np.ones(len(self.time))*Pos_ini.geocentric[0]
        y = np.ones(len(self.time))*Pos_ini.geocentric[1]
        z = np.ones(len(self.time))*Pos_ini.geocentric[2]
        
        return list([x,y,z])
    
  