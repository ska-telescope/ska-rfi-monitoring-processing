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
from scipy import signal

ms = 1e-3
us = 1e-6
MHz = 1e6
km = 1e3
minute = 60
hr = 60*minute
km_h = km/hr


#%% Cosmic background
'''
11/06/2019
Federico Di Vruno
This function generates time domain noise for low frequencies (<350 MHz).

Inputs:

Outputs:
    
'''    

def Cosmic_background(time_vect,fcentre,integ_time):
    f = fcentre/1e6
    #Tsky_ITU = 2.7+200*(408/f)**2.75  # per ITU-R P.372-7 page 19.
    
    Tsky = 35+26*(408/f)**2.75  # per Signal chain document for LOW.
    
    
    #1000 hs:
    B = 5400
    hs = 1000

    Tsky_ave = Tsky/np.sqrt(B*hs*3600)
    
    print('Ave Tsky @%.2f MHz = %.1f'%(f,Tsky))
    
    
    k = 1.38e-23
    B = 5400
    P_ave = 10*np.log10(k*Tsky_ave*B*1e3)
    print('1000 hs ave Tsky_ITU Pow @%.2f MHz = %.1f'%(f,P_ave-10))
    
    

#%% Signal
'''
11/06/2019
Federico Di Vruno
This function generates N periods of any signal selected, at the moment only ADS-B and DME.
Inputs:
        Name: signal name, at the moment only DME and ADS-B
        SamplingFreq: in Hz
        fcentre: centre frequency in Hz
        Ampl_change: is set to 1 randomize the amplitude from 0.1 to 1.1 with rectangular probability distribution
        duration: length of the signal generated in seconds.
        rand_phase: sets a random phase to the carrier signal with rect. PDF
        rand_displace: shifts the position of the pulse envelope by a random quantity within [0:period/2]

output:
        t_aux: time vector 
        S: voltage signal scaled to 1 (1.1 in the case of random amplitude)
        
'''
def Signal(Name, SamplingFreq, fcentre, duration, Ampl_change=0, rand_phase=0, rand_displace=0):

    
    if Name == 'DME':
        period = 1*ms 
        N_pulses = int(duration/period)+1
    if Name == 'ADS-B':
        period = 1*ms 
        N_pulses = int(duration/period)+1

    
    t_aux = []
    env_aux = []
    
    for i in range(N_pulses):
        
        if Name == 'DME':
            # signal between 1025 to 1150 MHZ
            period = 1*ms # periodicity of the signal
            offset = 12*us # separation between the pulses
            
            fc = fcentre
            Bw = 1*MHz
    
            t = np.linspace(-period/2,period/2,int(period*SamplingFreq))
    
            [k,env1] = Sig.gausspulse(t,fc,Bw/fc,retenv=1) #envolvent of the gauss pulse
    
            env2 = np.roll(env1,int(offset*SamplingFreq)) #envolvent shifted in time to generate the second pulse.
            
            t = np.linspace(0,period,int(period*SamplingFreq))
            env = env1 + env2
    
    
        if Name == 'ADS-B':
            # 
            Ton = 120*us
            fc = fcentre
            period = 1*ms # periodicity of the signal
            
            Trise = 0.1*us
            Tfall = 0.2*us
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
            t = np.linspace(0,period,int(period*SamplingFreq))
            
        try:
            t_aux = np.concatenate((t_aux,t_aux[-1]+t))
        except:
            t_aux = np.concatenate((t_aux,t))
                
        env_aux = np.concatenate((env_aux,env))


    if rand_phase:
        phase = np.random.rand(1)*2*np.pi
    else:
        phase = 0
        
    S = np.sin(2*np.pi*fc*t_aux+phase)*env_aux
    
    # displacement of the signal
    if rand_displace:
        N = int(fs/period)
        displace = int(np.random.rand(1)*N/2)
        S = np.roll(S,displace)
          
    if Ampl_change:
        A = np.random.rand(1)+0.1
        S = S*A
        
    return [t_aux[(t_aux<=duration)],S[(t_aux<=duration)]]



def Gauss_Noise(t,fs,stdev):
    # Define a white gausian noise to add to the RFI signal.
    # noise is defined in voltage.
    # Bandwidth of the noise is -fs/2 to fs/2
    N = len(t)
    Noise = stdev*np.random.randn(N)
    return Noise


'''

11/06/2019
Federico Di Vruno
This function generates N periods of any signal selected, at the moment only ADS-B and DME.

name: multiple_emitters

Inputs:
    Name: 
    N_emitters: 
    duration: time duration in sec
    fs: Sampling freq [Hz]

Output:
    t1: time vector
    S = signal
    '''

def multiple_emitters(Name,N_emitters,duration,fs):
    # This function generates a repetition of the basic signal N_repeat times,
    if Name == 'DME':
        fc = 1025*MHz + np.random.rand(N_emitters)*125*MHz # randomize the center frequency.
        
    if Name == 'ADS-B':
        fc = np.ones(N_emitters)*1030*MHz
        
    for i in range(N_emitters):
        t1,S1 = Signal(Name,fs,fc[i],duration,Ampl_change=1,rand_phase=1,rand_displace=0)
        if i==0: S = S1
        S += S1
    return t1[(t1<=duration)],S[(t1<=duration)]


def Tx_filter(S,fs,f1,f2):
#    fs = 800e6
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
    
    plt.figure()
    plt.plot( S, 'b', alpha=0.75)
    plt.plot( y, 'k')
    plt.legend(('noisy signal','filtfilt'), loc='best')
    
    [freq,V_orig] = RFI.fft_calc(S,fs, 1, 0)
    [freq,V_filt] = RFI.fft_calc(y,fs, 1, 0)
    
    plt.figure()
    plt.plot(freq,20*np.log10(V_orig),'r',freq,20*np.log10(V_filt),'b')

    return y
  
    
#%% DME
#
#[t1,S1] = Signal('DME',5000*MHz,1025*MHz)
##[t1,S2 ]= Signal('DME',5000*MHz)
##[t1,S3 ]= Signal('DME',5000*MHz)
#
#plt.figure()
#plt.plot(t1*1e6,S1)
##plt.plot(t1,S2)
##plt.plot(t1,S3)
#

#%%
fs = 3000*MHz
duration= 5*ms


#[t1,S_DME] = multiple_emitters('DME',20,duration,fs)
#plt.figure()
#plt.plot(t1,S_DME)


[t1,S_ADSB] = multiple_emitters('ADS-B',5,duration,fs)
plt.figure()
plt.plot(t1,S_ADSB)


#S = S_DME + S_ADSB
S = S_ADSB

S_filtered = Tx_filter(S,fs,(1030-20)*MHz,(1030+20)*MHz)

f_fft,P_fft = RFI.fft_calc(S,fs,1,0)
f_fft,P_fft_filt = RFI.fft_calc(S_filtered,fs,1,0)


#%% 
plt.figure()
plt.plot(f_fft,10*np.log10(P_fft_filt*1e3))
fo = 1030*MHz
f_ADSB_mask = np.array([-100,-78,-78,-23,-23,-7 ,-7,-1.3,-1.3,1.3,1.3,7 ,7  ,23 ,23 ,78 ,78 ,100])*MHz+fo
Amp = np.max(10*np.log10(P_fft*1e3))
Amp_mask = np.array([   -60 ,-60,-40,-40,-20,-20,-3,-3  ,0   ,0  ,-3 ,-3,-20,-20,-40,-40,-60,-60]) + Amp
plt.plot(f_ADSB_mask/1e6,Amp_mask,'r')




#%% Time occupied by the signal

count = np.sum((abs(S)>0.01))*100/len(S)


#%%
#def RFI_use_case_generator(signal_type,fs,duration,precision = 'Float64'):
    

