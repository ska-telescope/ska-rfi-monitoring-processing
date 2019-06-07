# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:46:31 2019

Signal generator for RFI use cases for SKA.


@author: f.divruno
"""
import scipy.signal as Sig
import numpy as np
import matplotlib.pyplot as plt
import RFI_general_functions as RFI

ms = 1e-3
us = 1e-6
MHz = 1e6
km = 1e3
minute = 60
hr = 60*minute
km_h = km/hr

#%%
def Signal(Name,SamplingFreq,rand_freq=0,rand_phase=0, rand_displace=0):
# this function generates one period of any signal.
    
    if Name == 'DME':
        # signal between 1025 to 1150 MHZ
        period = 1*ms # periodicity of the signal
        offset = 100*us # separation between the pulses
        
        fc = 1025*MHz
        if rand_freq: fc = fc + np.random.rand(1)*125 # randomize the center frequency.

        Bw = 1*MHz
        t = np.linspace(-period/2,period/2,period*SamplingFreq)

        [k,env1] = Sig.gausspulse(t,fc,Bw/fc,retenv=1) #envolvent of the gauss pulse
        #env1 = np.roll(env1,int(-period/2*SamplingFreq))
        env2 = np.roll(env1,int(offset*SamplingFreq)) #envolvent shifted in time to generate the second pulse.
        
        t = np.linspace(0,period,period*SamplingFreq)
        env = env1 + env2


    if Name == 'ADS-B':
        # 
        Ton = 120*us
        fc = 1090*MHz #1090*MHz
        period = 1*ms # periodicity of the signal
        
        Trise = 0.1*us
        Tfall = 0.1*us
        Tpulse = 0.5*us-Trise-Tfall
        pulse_env = np.concatenate((np.linspace(0,1,int(Trise*SamplingFreq)) \
                                    ,np.ones(int(Tpulse*SamplingFreq)) \
                                    ,np.linspace(1,0,int(Tfall*SamplingFreq)) \
                                    ,np.zeros(int((1*us-Trise-Tfall-Tpulse)*SamplingFreq))))
        env = []
        for i in range(120):
            pulse_shift = (np.round(np.random.random(1))).astype(int)
            env = np.concatenate((env,np.roll(pulse_env,int(pulse_shift*0.5*us*SamplingFreq))))
        
        env = np.concatenate((env,np.zeros(int((period-Ton)*SamplingFreq))))
        t = np.linspace(0,period,period*SamplingFreq)
        

    if rand_phase:
        phase = np.random.rand(1)*2*np.pi
    else:
        phase = 0
    sine = np.sin(2*np.pi*fc*t+phase)
    S = sine*env
    
    # displacement of the pulses inside the time vector.
    if rand_displace:
        N = len(S)
        displace = int(np.random.rand(1)*N/2)
        S = np.roll(S,displace)
            
                
    return [t,S]


  
    
#%% DME

[t1,S1] = Signal('DME',5000*MHz)
#[t1,S2 ]= Signal('DME',5000*MHz)
#[t1,S3 ]= Signal('DME',5000*MHz)

plt.figure()
plt.plot(t1,S1)
#plt.plot(t1,S2)
#plt.plot(t1,S3)


#[t1,S_DME] = gen_signal_train('DME',5000*MHz,5)
#plt.figure()
#plt.plot(t1,S_DME)



#%%
def RFI_use_case_generator(signal_type,fs,duration,precision = 'Float64'):
    

