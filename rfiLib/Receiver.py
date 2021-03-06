# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 23:50:54 2019

@author: f.divruno
@revised: G. Hovey; added default parameters to initialize class
@revised: G. Hovey; somehow a copy with antenna_gain unindented got into the repository causing class functions to disappear. Fix it.
@revised: G. Hovey; changed height to 1000 u*m from zero (this needs to be a parameter)
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as Coord
import astropy.units as u
from rfiLib.adcGain import adcGain
from rfiLib.siggen_aph import WhiteNoiseSignal, band_limit

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
    def __init__(self,Name='',
                 Position=dict(Latitude = -30.71329*u.deg, \
                               Longitude = 21.449412*u.deg),\
                               Pointing=dict(Elev=90*u.deg,Azimuth=0*u.deg),\
                               duration=.5*ms,SampleRate=4.*GHz,antSampleRate=3.96*GHz,\
                               band = 'B2'):
        self.band = band
        self.height = 1000 *u.m
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
        self.antSampleRate = antSampleRate
        self.Duration = duration
        self.sky_source_rx = []
        self.Rx_signal = []
        self.time = []
        self.sky_source_rx_foffset = []
        self.Rx_signal_foffset = []
        self.time_foffset = []

        #self.Pointing = self.Pointing()
        
        
        
    def antenna_gain(angle):
        
        N = np.size(angle)
        G = np.zeros(N)
        delay_samples = np.zeros(N)
        angle = np.abs(angle)
        if N >1:
            for i in range(N):
                if angle[i]<0.1:
                    G[i] = 60
                    delay_samples[i] = 0
            
                elif (angle[i]>0.1) & (angle[i]<3.5):
                    G[i] = 20-40*np.log10(angle[i]) # dB
                    delay_samples[i] = np.random.rand(1)*10 # maximum random delay because of the sidelobes = 10/fs = 2.5 ns en 4Gsps
                        
                else:
                    G[i]=0
                    delay_samples[i] = np.random.rand(1)*10
        else:
            if angle<0.1:
                G = 60
                delay_samples = 0
        
            elif (angle>0.1) & (angle<3.5):
                G = 32-25*np.log10(angle) # dB
                delay_samples = np.random.rand(1)*10 # maximum random delay because of the sidelobes = 10/fs = 2.5 ns en 4Gsps
                    
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
    
    def plot_signal(self,tap,signal,mode):
        '''
           tap: antIn, adcIn, adcOut 
           mode: absVolt, volt, powerLin, powerLog
           signal: RFI,Sky
           Plots the time domain signal in power or voltage 
           
        '''
        if tap == 'antIn':
            if signal == 'RFI':
                sig = np.copy(self.Rx_signal)
            else:
                sig = np.copy(self.sky_source_rx)
        elif tap == 'adcIn':
            if signal == 'RFI':
                sig = np.copy(self.ADC_input_rx)
            else:
                sig = np.copy(self.ADC_input_sky)
        elif tap == 'adcOut':
            if signal == 'RFI':
                sig = np.copy(self.ADC_output_rx)
            else:
                sig = np.copy(self.ADC_output_sky)
        else:
            print('The tap selected is non existant: \nChoose: antIn, adcIn, adcOut')
            raise SystemExit    
    
        time = self.time        

               
        if mode == 'absVolt':
            S = np.abs(sig)         #Instantaneous voltage in absolute V
            unit = 'magnitude (V)'
        
        elif mode == 'volt':
            S = sig                 #Instantaneous voltage in V
            unit = 'V'
        
        elif mode =='powerLin':
            S = np.abs(sig)**2/50   # Instantaneous power in W
            unit = 'Watts'
            
        elif mode =='powerLog':     # Instantaneous power in dBm
            S = 10*np.log10(np.abs(sig)**2/50*1e3) 
            unit = 'dBm'
            
        elif mode == 'voltLog':     # Instantaneous Voltage in dBV
            S = 20*np.log10(np.abs(sig))
            unit = 'dBV'


        if tap == 'adcOut':
            S = sig         #in ADC out the value is ADC counts
            unit = 'ADC counts'
        
        plt.figure()
        plt.plot(time/us,S)
        plt.title('%s signal in %s,antenna %s'%(signal,tap,self.Name))
        plt.ylabel(unit)
        plt.xlabel('us')

    def plot_spectrum(self,tap,signal, mode='power'):
        '''
            tap: antIn, adcIn, adcOut 
            mode: power, voltage
            signal: RFI , Sky
        '''
        if tap == 'antIn':
            if signal == 'RFI':
                sig = np.copy(self.Rx_signal)
            else:
                sig = np.copy(self.sky_source_rx)
        elif tap == 'adcIn':
            if signal == 'RFI':
                sig = np.copy(self.ADC_input_rx)
            else:
                sig = np.copy(self.ADC_input_sky)
        elif tap == 'adcOut':
            if signal == 'RFI':
                sig = np.copy(self.ADC_output_rx)
            else:
                sig = np.copy(self.ADC_output_sky)
        else:
            print('The tap selected is non existant: \nChoose: antIn, adcIn, adcOut')
            raise SystemExit    


        V = sig
        fs = self.SampleRate
        N = len(V)
        if tap == 'adcOut':
            S = 10*np.log10(np.abs(np.fft.fft(V)/N))
            ylabel = 'Counts [dBcounts]'
        else:    
            if mode == 'power':
                S = 10*np.log10(np.abs(np.fft.fft(V)/N)**2/50*1e3)
                ylabel = 'Power [dBm]'
            else:
                S = 10*np.log10(np.abs(np.fft.fft(V)/N)**2)
                ylabel = 'Voltage [dBV]'

        freq = np.fft.fftshift(np.fft.fftfreq(N, d=1/fs))
        S = np.fft.fftshift(S) 
        
        plt.figure()
        plt.plot(freq/MHz,S)
        plt.title('%s signal in %s,antenna %s'%(signal,tap,self.Name))
        plt.ylabel(ylabel)
        plt.xlabel('MHz')
        
    

    def Apply_antSampleRate(self): 
        """
            Re-sampling by linear interpolation (only down-sampling). Both time and amplitudes are re-sampled
            without scale change, so ENERGY IS CONSERVED but POWER IS NOT.
            @param f_s: sampling rate for "signal"
            @param f_s_new: desired sampling rate, must be <= f_s.
            @return: (t_new, s_new) re-sampled copies of time and signal series.
             by Adriaan Peens-Hugh
             modified by FDV 30/07/19 
        """
        
        f_s_new = self.antSampleRate
        
        t_new = np.linspace(1/f_s_new, self.time[-1],self.time[-1]*f_s_new )

        self.ADC_output_rx_foffset = np.interp(t_new, self.time, self.ADC_output_rx)
        
        self.ADC_output_sky_foffset = np.interp(t_new, self.time, self.ADC_output_sky)
        self.time_foffset = t_new

    
    def Noise_scale(self, s,scaling,nBits,Gain = 0, Vfs=1):
        """
            returns a gain in dB to scale the noise signal, so that 1 digitizer level is
            E = 0.335*Vrms in the case of Correlator_optimized.
            In the case of Linearity_optimized scales the peak signal to fit in the ADC range.
           
        """
        if scaling == 'Correlator_opimized':
            Vrms = np.sqrt(np.sum(s**2)/len(s))
            PdBm = 10*np.log10(Vrms**2/50*1e3)
            Gain = adcGain(PdBm,nBits,Vfs) #dB gain to scale the input voltage.
            
            
        elif scaling == 'Linearity_opimized':
            Vrms = np.sqrt(np.sum(s**2)/len(s))
            Vpk = Vrms*2**0.5
            Gain = 20*np.log10(Vfs/Vpk)
    
        elif scaling == 'Defined_Gain':
            Gain = Gain
            
        self.Analog_Gain_dB = Gain    
        self.Analog_Gain_lin = 10**(Gain/10)

        return Gain
    
    def quantize(self,signal, nbits, Vfs=1):
        """ 
            Quantizes the input signal to discrete, equi-spaced levels. Even though the output are
            discrete values representable by integers from 0 to full scale
            @param level_width: typically N*std(signal) with N = 0.335 for 4b [EVLA memo 88].
            @return: the integer quantized version of the original "analogue" signal.
            by Adriaan Peens-Hugh
        """
        bins = np.arange(-Vfs,Vfs,Vfs/(2**(nbits-1)))
        
        signal[signal<bins[0]] = bins[0] # (-inf,min]   #Clippiong to the minimum level    
        signal[signal>bins[-2]] = bins[-1] # (max,+inf)  # clipping to the maximum level 
    
    
        qsignal = np.digitize(signal,bins,True) - 2**(nbits-1)
        
        
        # Convert to integers from 0 to fullscale
        return qsignal.astype('int')    
    
    def Apply_ADC(self,nBits):
        '''
            Applies the quantizer to the signal (already amplified)
        '''
        _ADC_input_rx = np.copy(self.ADC_input_rx)
        _ADC_input_sky = np.copy(self.ADC_input_sky)
        self.ADC_output_rx = self.quantize(_ADC_input_rx, nBits )
        self.ADC_output_sky = self.quantize(_ADC_input_sky, nBits)
        
        
    def Apply_analog_chain(self,Band,scaling,atten=0,f_offset=0):
        '''
            Applies the zero order model of the analog chain of DISH
        '''
        s_sky = np.copy(self.sky_source_rx) #signal with only sky source
        s_rx = np.copy(self.Rx_signal) # signal with sky source and RFI
        t_s = np.copy(self.time)
        f_s = np.copy(self.SampleRate)
        
        if Band == 'B1':
            print('Applying Band 1')
            """Zero'th order model for MID B1 EM+analogue section """
            # for the signal only case
            s = s_sky
            # Approx. T_ant = sky + atmosphere + spillover @ ZA=0 (40@350MHz, 10@1050MHz) [SKA-TEL-DSH-0000082 rev 1, fig 20]
            tAntNoise = WhiteNoiseSignal(t_s, Teq=25)
            # Approx. T_rx (18.5@1050MHz) [SKA-TEL-DSH-0000082, rev 1, table 2]
            tRxNoise = WhiteNoiseSignal(t_s, Teq=18.5)
            
            s += tAntNoise + tRxNoise
            # Filtering
            s = band_limit(s, f_s, (350*MHz,1050*MHz), (300*MHz,1200*MHz), 0.5, 40, ftype="ellip") # [SKA-TEL-DSH-00000021 rev 2, fig 14]
            
            # Attenuation scaling, according to the scaling strategy
            Gain = self.Noise_scale(s,scaling,12)
            self.Atten =  Gain-56-32.5 # calculates the real value needed to set the attenuators.
            
            # Gain
            s *= 10**((56+32.5+self.Atten)/20.)
            s_sky = np.copy(s)



            # for the signal + RFI case
            s = s_rx
            
            # Approx. T_ant = sky + atmosphere + spillover @ ZA=0 (40@350MHz, 10@1050MHz) [SKA-TEL-DSH-0000082 rev 1, fig 20]
            s += tAntNoise + tRxNoise
            # Filtering
            s = band_limit(s, f_s, (350*MHz,1050*MHz), (300*MHz,1200*MHz), 0.5, 40, ftype="ellip") # [SKA-TEL-DSH-00000021 rev 2, fig 14]
            # Gain
            s *= 10**((56+32.5+self.Atten)/20.)
            s_rx = np.copy(s)




        if Band == 'B2':
            print('Applying Band 2')
            """Zero'th order model for MID B2 EM+analogue section """
            # for the signal only case
            s = np.copy(s_sky)
            # Approx. T_ant = sky + atmosphere + spillover ZA=0 (8.2@950MHz, 5.8@1760MHz) [SKA-TEL-DSH-0000111, rev 1, fig 37]
            tAntNoise = WhiteNoiseSignal(t_s, Teq=7)
            # Approx. T_rx (10.6@1650MHz, 12.1@3050MHz) SKA-TEL-DSH-0000111, rev 1
            tRxNoise = WhiteNoiseSignal(t_s, Teq=11)
            
            s += tAntNoise + tRxNoise
            # Filtering
            s = band_limit(s, f_s, (950*MHz,1760*MHz), (800*MHz,2000*MHz), 0.5, 50, ftype="cheby1") # [SKA-TEL-DSH-00000021 rev 2, fig 15]
             
            # Attenuation scaling, according to the scaling strategy
            Gain = self.Noise_scale(s,scaling,12)
            self.Atten =  Gain-56-32 # calculates the real value needed to set the attenuators.
            
            # Gain
            s *= 10**((56+32+self.Atten)/20.)
            s_sky = np.copy(s)                       
            
            # for the signal + RFI case            
            s = np.copy(s_rx)
            s += tAntNoise + tRxNoise

            # Filtering
            s = band_limit(s, f_s, (950*MHz,1760*MHz), (800*MHz,2000*MHz), 0.5, 50, ftype="cheby1") # [SKA-TEL-DSH-00000021 rev 2, fig 15]
             
            # Gain
            s *= 10**((56+32+self.Atten)/20.)
            s_rx = np.copy(s)


            
        if Band == 'B3':
            """Zero'th order model for MID B3 EM+analogue section """
            # Approx. T_ant = sky + atmosphere + spillover @ ZA=0 (8@1650MHz, 6.5@3050MHz) [SKA-TEL.DSH.SE-NRF-R-001 rev 2, p 202]
            s = s + WhiteNoiseSignal(t_s, Teq=7)
            # Approx. T_rx (10.6@1650MHz, 12.1@3050MHz)-2 for DS [SKA-TEL.DSH.SE-NRF-R-001 rev 2, table 8]
            s += WhiteNoiseSignal(t_s, Teq=10)
            # Filtering
            s = band_limit(s, f_s, (1650*MHz,3050*MHz), (1500*MHz,3250*MHz), 0.5, 50, ftype="ellip") # [SKA-TEL-DSH-00000021 rev 2, fig 16]
            
            # Attenuation scaling, according to the scaling strategy
            Gain = self.Noise_scale(s,scaling,12)
            self.Atten =  Gain+56+29 # calculates the real value needed to set the attenuators.
            
            # Gain
            s *= 10**((56+29-self.Atten)/20.)

        if Band == 'B4':
            """Zero'th order model for MID B4 EM+analogue section """
            # Approx. T_ant = sky + atmosphere + spillover @ ZA=0 (7.5+-1.5@2800MHz-5180MHz) [SKA-TEL.DSH.SE-NRF-R-001 rev 2, p 242]
            s = s + WhiteNoiseSignal(t_s, Teq=7.5)
            # Approx. T_rx (14.3@2800MHz, 16.7@5180MHz)-2 for DS [SKA-TEL.DSH.SE-NRF-R-001 rev 2, table 10]
            s += WhiteNoiseSignal(t_s, Teq=12)
            # Filtering
            s = band_limit(s, f_s, (2800*MHz,5180*MHz), (2500*MHz,5800*MHz), 0.5, 60, ftype="cheby1") # [318-000000-005, rev A, fig 6a]
            
            # Attenuation scaling, according to the scaling strategy
            Gain = self.Noise_scale(scaling,12)
            self.Atten =  Gain+56+18+6 # calculates the real value needed to set the attenuators.
            
            # Gain
            s *= 10**((56+18+6-self.Atten)/20.)


        if Band == 'B5a':
            """Zero'th order model for MID B5a EM+analogue section """
            # Approx. T_ant = sky + atmosphere + spillover @ ZA=0 (7+-1@4600MHz-8500MHz) [SKA-TEL-DSH-0000118, rev 1, fig 76]
            s = s + WhiteNoiseSignal(t_s, Teq=7)
            # Approx. T_rx (7.7+0.633/GHz*2GHz@65500MHz) - 2.5K for DS [SKA-TEL-DSH-0000118, rev 1, table 26]
            s += WhiteNoiseSignal(t_s, Teq=6.5)
            # Filtering
            s = band_limit(s, f_s, (4600*MHz,8500*MHz), (4000*MHz,9500*MHz), 0.5, 60, ftype="cheby1") # [SRxB45, fig 6b]
            
            # Attenuation scaling, according to the scaling strategy
            Gain = self.Noise_scale(s,scaling,12)
            self.Atten =  Gain+56+16+8 # calculates the real value needed to set the attenuators.
            
            # Gain
            s *= 10**((56+16+8-self.Atten)/20.)

                    
        if Band == 'B5b':            
            """Zero'th order model for MID B5b EM+analogue section """
            # Approx. T_ant = sky + atmosphere + spillover @ ZA=0 (8.7+-1@8300MHz-15400MHz) [SKA-TEL-DSH-0000118, rev 1, fig 76]
            s = s + WhiteNoiseSignal(t_s, Teq=8.7)
            # Approx. T_rx (10.6+0.633/GHz*3GHz@11410MHz) - 2.5K for DS [SKA-TEL-DSH-0000118, rev 1, table 27]
            s += WhiteNoiseSignal(t_s, Teq=10)
            # filtering
            s = band_limit(s, f_s, (8300*MHz,15400*MHz), (6500*MHz,19500*MHz), 0.5, 60, ftype="cheby1") # [SRxB45, fig 6c]
                        
            # Attenuation scaling, according to the scaling strategy
            Gain = self.Noise_scale(s,scaling,12)
            self.Atten =  Gain+56+14+4 # calculates the real value needed to set the attenuators.
            
            # Gain
            s *= 10**((56+14+4-self.Atten)/20.)
            
        self.ADC_input_rx = np.copy(s_rx)    #Sky ignal plus RFI
        self.ADC_input_sky = np.copy(s_sky)   # only sky signal and white noise from rcvr.
        