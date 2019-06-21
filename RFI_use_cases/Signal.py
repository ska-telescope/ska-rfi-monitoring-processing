# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 23:47:19 2019

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
from scipy import signal

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
class Signal():

    # generate a signal object.
    def __init__(self,Name,Duration,SamplingRate,CenterFreq,Power,Pol,random_seed=[]):
        self.Name = Name
        self.Duration = Duration
        self.SamplingRate = SamplingRate
        self.CenterFreq = CenterFreq
        self.Seed = random_seed #Used in the case of wanting to repeat the "random" signal (not implemented yet)
        self.Power = Power
        self.Polarization = Pol
        
        self.time,self.data = self.Message()
        
        self.data = self.Scale_power()

    
 # --------------------------------------        
    def Symbol(self, ini_phase, displace = 0):
        
        if self.Name == 'DME':
            period = 1*ms # periodicity of the signal
            offset = 12*us
            SamplingFreq = self.SamplingRate
            fc = self.CenterFreq
            Bw = 0.5*MHz #0.5% of center frequency
            t = np.linspace(-period/2,period/2,period*SamplingFreq)
    
            [k,env1] = Sig.gausspulse(t,fc,(Bw/fc),retenv=1) #envelope of the gauss pulse
            env2 = np.roll(env1,int(offset*SamplingFreq)) #envelope shifted in time to generate the second pulse.
            t = np.linspace(0,period,period*SamplingFreq)
            env = env1 + env2        
        
            sine = np.sin(2*np.pi*fc*t+ini_phase)
            S = sine*env
            
            # displacement of the pulses inside the time vector.
            if displace>0:
                #np.random.seed(self.Seed) #loads the seed in case of wanting to repeat the signal
                S = np.roll(S,displace)
            return [t,S]
        

        if self.Name == 'ADS-B':
             
            Ton = 120*us
            fc = self.CenterFreq #1090*MHz
            period = 1*ms # periodicity of the signal
            SamplingFreq = self.SamplingRate
            Trise = 0.1*us
            Tfall = 0.1*us
            Tpulse = 0.5*us-Trise-Tfall
            pulse_env = np.concatenate((np.linspace(0,1,int(Trise*SamplingFreq)) \
                                        ,np.ones(int(Tpulse*SamplingFreq)) \
                                        ,np.linspace(1,0,int(Tfall*SamplingFreq)) \
                                        ,np.zeros(int((1*us-Trise-Tfall-Tpulse)*SamplingFreq))))
            N = len(pulse_env)
            env = np.zeros(int(period*SamplingFreq))
            for i in range(120):
                try:
                    np.random.seed(self.Seed*i)
                except:
                    np.random.seed()
                pulse_shift = (np.round(np.random.random(1))).astype(int)
                env[i*N:(i+1)*N] = np.roll(pulse_env,int(pulse_shift*period/2*SamplingFreq))
            
            t = np.linspace(0,period,period*SamplingFreq)
            
            sine = np.sin(2*np.pi*fc*t+ini_phase)
            S = sine*env
        
            # displacement of the pulses inside the time vector.
            if displace>0:
                S = np.roll(S,displace)

            return [t,S]    
# --------------------------------------        
    def Message(self):
        L = self.Duration
        fs = self.SamplingRate
        
        
        if self.Name == 'DME':
            SymLen = 1e-3 #lenght of the symbol
            samples_sym = int((SymLen*fs)) #number of samples in a symbol
            N_Symbols = int((L/SymLen)) #number of symbles + 1 because of the rounding.
            samples_tot = int((N_Symbols*samples_sym)) # total number of samples
            #np.random.seed(self.Seed) #loads the seed in case of wanting to repeat the same signal.
            ini_phase = np.random.rand(1)*2*np.pi #random phase for the signal generator
            S = np.zeros(samples_tot)
            t = np.zeros(samples_tot)

            N = samples_sym
            displace = int(np.random.rand(1)*(N*0.25))
            for i in range(N_Symbols):
                t_aux,S[i*samples_sym:(i+1)*samples_sym] = self.Symbol(ini_phase, displace)        
                if i>0:
                    t_aux += t_aux[1]

                t[i*samples_sym:(i+1)*samples_sym] = t_aux+t_aux[-1]*i
            return t[t<=self.Duration] , S[t<=self.Duration]
            
        if self.Name == 'ADS-B':
            SymLen = 1*ms #time duration of a symbol >120us
            samples_sym = int(SymLen*fs) #number of samples in a symbol
            N_Symbols = int(L/SymLen)+1 #number of symbles + 1 because of the rounding.
            samples_tot = int(round(N_Symbols*SymLen*fs)) # total number of samples
            #np.random.seed(self.Seed) #loads the seed in case of wanting to repeat the same signal.
            ini_phase = np.random.rand(1)*2*np.pi #random phase for the signal generator
            N = samples_sym
            displace = int(np.random.rand(1)*(N*0.25))
            S = np.zeros(samples_tot)
            t = np.zeros(samples_tot)
            for i in range(N_Symbols):
                t_aux,S[i*samples_sym:(i+1)*samples_sym] = self.Symbol(ini_phase, displace)        
                if i>0:
                    t_aux += t_aux[1]

                t[i*samples_sym:(i+1)*samples_sym] = t_aux+t_aux[-1]*i
            S = self.Tx_filter(S,self.CenterFreq*0.99,self.CenterFreq*1.01)
            return t[t<=self.Duration] , S[t<=self.Duration]

    '''----------------------------------
    
        ----------------------------------
    '''
    
    def Scale_power(self):
        Pow = self.Power # power in dBm
        scale = np.sqrt(10**(Pow/10)*1e-3*50*2)
        return self.data*scale


    def Tx_filter(self,S,f1,f2,plot=0):
        fs = self.SamplingRate
        fn = fs/2
        f1 = f1/fn
        f2 = f2/fn
    
        b, a = signal.butter(3, [f1,f2],btype='bandpass',output='ba')
        
    #    zi = signal.lfilter_zi(b, a)
    #    z, _ = signal.lfilter(b, a, S, zi=zi*S[0])
    #    #Apply the filter again, to have a result filtered at an order the same as filtfilt:
    #    
    #    z2, _ = signal.lfilter(b, a, z, zi=zi*z[0])
    #    #Use filtfilt to apply the filter:
        y = signal.filtfilt(b, a, S)
        
        if plot:
            plt.figure()
            plt.plot( S, 'b', alpha=0.75)
            plt.plot( y, 'k')
            plt.legend(('noisy signal','filtfilt'), loc='best')
            
            [freq,V_orig] = RFI.fft_calc(S,fs, 1, 0)
            [freq,V_filt] = RFI.fft_calc(y,fs, 1, 0)
            
            plt.figure()
            plt.plot(freq,20*np.log10(V_orig),'r',freq,20*np.log10(V_filt),'b')
    
        return y