# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 23:50:54 2019

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
        
class Receiver():
    def __init__(self,Name,Position,Pointing,duration,SampleRate):
        self.height = 0
        self.Name = Name
        self.lat = Position['Latitude']
        self.lon = Position['Longitude']
        self.Coordinates = Coord.EarthLocation.from_geodetic(self.lon,self.lat,self.height) # position in XYZ
        Posx = self.Coordinates.geocentric[0].to(u.km)
        Posy = self.Coordinates.geocentric[1].to(u.km)
        Posz = self.Coordinates.geocentric[2].to(u.km)
        self.Pos = [Posx,Posy,Posz]
        self.Pointing = Pointing
        self.SampleRate = SampleRate
        #self.Pointing = self.Pointing()
        
        
        
    def antenna_gain(angle):
        if angle<1:
            G = 60
            delay_samples = 0

        elif (angle>1) & (angle<48):
            G = 32-25*np.log10(angle)
            delay_samples = np.random.rand(1)*10 # maximum random delay because of the sidelobes = 10/fs = 2.5 ns en 4Gsps
            return G            
        else:
            G=0
            delay_samples = np.random.rand(1)*10
        
        return delay_samples,G # gain in the beam    
    
    def plot_gain(self):
        for angle in range(180):
            delay,G = self.antenna_gain(angle)
        plt.figure()
        plt.plot(G)
        plt.title('Gain')
        plt.xlabel('')

    
    def Pointing():
        elev = 90
        azimuth = 0
        return elev,azimuth
    
    def plot_signal(self,mode,signal):
        '''
            mode: abs,abslog,lin
            signal: RFI,Sky
        '''
        if signal == 'RFI':
            sig = self.Rx_signal
            time = self.time
        else:
            sig = self.sky_source_rx
            time = np.linspace(0,self.SampleRate*len(sig),len(sig))
               
        if mode == 'abs':
            S = np.abs(sig)
        
        if mode == 'abslog':
            S = 20*np.log10(np.abs(sig))

        if mode == 'lin':
            S = 20*np.log10(np.abs(sig))

        
        plt.figure()
        plt.plot(time/us,S)
        plt.title('Received signal in telescope '+self.Name)
        plt.ylabel(mode)
        plt.xlabel('us')

    def plot_spectrum(self,mode,signal):
        '''
            mode: 
            signal: RFI , Sky
        '''
        if signal == 'RFI':
            sig = self.Rx_signal
        else:
            sig = self.sky_source_rx

        V = sig
        fs = self.SampleRate
        N = len(V)
        S = 10*np.log10(np.abs(np.fft.fft(V)))
        freq = (np.fft.fftfreq(N, d=1/fs))
        
        
        plt.figure()
        plt.plot(freq/MHz,S)
        plt.title('Received %s signal in telescope '%(signal)+self.Name)
        plt.ylabel('dB')
        plt.xlabel('MHz')
        
    

    def DownSample(signal, f_s, f_s_new):
        """
            Re-sampling by linear interpolation (only down-sampling). Both time and amplitudes are re-sampled
            without scale change, so ENERGY IS CONSERVED but POWER IS NOT.
            @param f_s: sampling rate for "signal"
            @param f_s_new: desired sampling rate, must be <= f_s.
            @return: (t_new, s_new) re-sampled copies of time and signal series.
             by Adriaan Peens-Hugh
        """
        assert (f_s >= f_s_new), "sample() does not support up-sampling!"
        t = np.arange(0, len(signal), step=1)*1/f_s
        t_new = np.arange(0, len(signal), step=f_s/f_s_new)*1/f_s
        signal = np.interp(t_new, t, signal)
        return t_new, signal
    
    def Noise_scale(self,scaling,Gain):
        """
            returns a gain in dB to scale the noise signal, so that 1 digitizer level is 
            
        """
        if scaling == 'Correlator_opimized':
            ####
            #### Insert Gary code here
            ####
            ####
            Gain =10
            
        elif scaling == 'Linearity_opimized':
            ####
            #### Insert Gary code here
            ####
            ####
            Gain = 10
    
        elif scaling == 'Defined_Gain':
            
        self.Analog_Gain_dB = Gain    
        self.Analog_Gain_lin = 10**(Gain/10)
    
    
    def Quantize(signal, nbits, Vfs):
        """ 
            Quantizes the input signal to discrete, equi-spaced levels. Even though the output are
            discrete values representable by integers from 0 to full scale
            @param level_width: typically N*std(signal) with N = 0.335 for 4b [EVLA memo 88].
            @return: the integer quantized version of the original "analogue" signal.
            by Adriaan Peens-Hugh
        """
        qsignal = np.copy(signal)
        yrange = np.arange(-2**(nbits-1),2**(nbits-1),1)
        yrange = yrange - np.mean(yrange) # Central value of threshold levels
        xthresh = yrange + 0.5 # Upper value of level
        qsignal[qsignal<=xthresh[0]] = yrange[0] # (-inf,min]   #Clippiong to the minimum level
        
        for l,u,y in zip(xthresh[:-2],xthresh[1:-1],yrange[1:-1]):
            qsignal[np.logical_and(qsignal>l,qsignal<=u)] = y # (-0.5,0.5]
            
        qsignal[qsignal>xthresh[-2]] = yrange[-1] # (max,+inf)  # clipping to the maximum level 
        # Convert to integers from 0 to fullscale
        return qsignal.astype('int')
    
