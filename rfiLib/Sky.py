# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 23:49:31 2019

@author: f.divruno
"""
import numpy as np
import astropy.coordinates as Coord
import astropy.constants as const

import rfiLib.General as General

import pyproj as pyproj
import astropy.units as u
from rfiLib.Pointing_DISH import  Pointing_to_ECEF
from rfiLib.siggen_aph import WhiteNoiseSignal

ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

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
    def __init__(self,Name,srcPointing,posRx, SamplingRate,Duration,Temperature,random_seed = None):
        """            
            srcPointing: dict(elev=...*u.deg,az=...*u.deg) 
            SamplingRate:
            Duration:
            Temperature:
            random_seed:
                
        """
        self.Duration = Duration
        self.SamplingRate = SamplingRate
        self.Temp = Temperature
        self.random_seed = random_seed
        self.BW = SamplingRate/2
        self.Power = k_bolt*Temperature*self.BW
        self.Ampl = np.sqrt(self.Power/50)
        self.Name = Name
               
        self.elev = srcPointing['elev']
        self.az = srcPointing['az']
  
        
        lonRx,latRx,altRx = pyproj.transform(ecef,lla,posRx[0].value*1e3,posRx[1].value*1e3,posRx[2].value*1e3) #plot all sats          
        
        self.srcECEF = Pointing_to_ECEF(self.elev,self.az,latRx*u.deg,lonRx*u.deg)*np.linalg.norm([posRx[0].value,posRx[1].value,posRx[2].value]) #ECEF versor pointing to the source
        
#        self.Pos = Coord.EarthLocation.from_geodetic(srcPosition['lon'],srcPosition['lat'])
        
        self.time, self.data = self.sky_signal()
 
        
        
    def sky_signal(self):
        # Original White noise signal to be summed to every antenna    
        # Sky is generated with an extra time to be able to delay the signal
        # at a maximum angle of 30 deg for a baseline of 80km
#        max_baseline = 80*km
        max_angle = 30*np.pi/180
#        max_delay = max_baseline/const.c*np.sin(max_angle)
        t_max = self.Duration
        
        #number of samples
        N = int(t_max*self.SamplingRate)
        time = np.linspace(1/self.SamplingRate,t_max,N)
        #generate the random signal
        sky = WhiteNoiseSignal(time,Teq=self.Temp,rand_seed= self.random_seed)
#        try:
#            sky = WhiteNoiseSignal(time,self.Temp,self.random_seed)
#        except:
#            sky = WhiteNoiseSignal(time,self.Temp,self.random_seed)
        print('Noise signal PSD %f'%(10*np.log10(np.std(sky)**2/50)+30-10*np.log10(self.SamplingRate/2)))
        return time,sky
        
    
    # Delay of the star signal depending on the baseline and the angle of the source    
    def rx_sky(self,SamplingRate,duration,posRx):
        
#        posSrc = General.Coord_to_nparray(self.Pos.geocentric)
        posSrc = self.srcECEF
        posRx = np.array(posRx)/km

   
        #delay:
        #the delay changes with the angle of incidence and the height of the transmiter.
        lonTx,latTx,altTx = pyproj.transform(ecef,lla,posSrc[0]*1e3,posSrc[1]*1e3,posSrc[2]*1e3) #plot all sats    
        lonRx,latRx,altRx = pyproj.transform(ecef,lla,posRx[0]*1e3,posRx[1]*1e3,posRx[2]*1e3) #plot all sats    
        
        
        #calculate the point on the earth where the phase reference is zero:
        posTxPhaseRef = np.array(pyproj.transform(lla,ecef,lonTx,latTx,altRx))/km
        distRxPhaseRef = np.linalg.norm(posTxPhaseRef-posRx)
        
        #calculate the angle of incidence
        dRxPhaseRef = posSrc-posRx
        theta0 = (np.arccos(np.dot(posSrc,posRx)/np.linalg.norm(posSrc)/np.linalg.norm(posRx))*180/np.pi)*u.deg
        theta1 = (np.arccos(np.dot(dRxPhaseRef,posRx)/np.linalg.norm(dRxPhaseRef)/np.linalg.norm(posRx))*180/np.pi)*u.deg
        alpha = 180*u.deg- theta0 - theta1
        
        #calculate the distance with respect to the phase centre
        d = np.linalg.norm(distRxPhaseRef)*np.cos(alpha)*u.km
        time_delay = ((d/const.c).to(u.s)).value
            

        delaySamples = int(time_delay*SamplingRate)

        
#        sky_delayed = self.data[N_delay:(N_delay+N)]
        
        return delaySamples


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
        
