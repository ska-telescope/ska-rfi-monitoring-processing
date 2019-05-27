# -*- coding: utf-8 -*-
"""
Created on Fri May 24 14:18:19 2019

@author: F.Divruno
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


def Signal(Name,SamplingFreq,rand_phase=0, rand_displace=0):
    
    if Name == 'DME':
        period = 1*ms # periodicity of the signal
        offset = 100*us
        fc = 100*MHz
        Bw = 1*MHz
        t = np.linspace(-period/2,period/2,period*SamplingFreq)

        [k,env1] = Sig.gausspulse(t,fc,Bw/fc,retenv=1) #envolvent of the gauss pulse
        env2 = np.roll(env1,int(offset*SamplingFreq)) #envolvent shifted in time to generate the second pulse.
        env = env1 + env2

    if Name == 'ADS-B':
        # 
        Ton = 120*us
        fc = 100*MHz #1090*MHz
        period = 1*ms # periodicity of the signal
        
        Trise = 0.05*us
        Tfall = 0.05*us
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

def Gauss_Noise(t,fs,sig,mu):
    # Define a white gausian noise to add to the RFI signal.
    # noise is defined in voltage.
    # Bandwidth of the noise is -fs/2 to fs/2
    N = len(t)
    Noise = sig*np.random.randn(N) + mu
    return Noise

    

def gen_signal_train(Name,fs,N_repeat):
    # This function generates a repetition of the basic signal N_repeat times,
    
    tout = []
    Sout = []

    for i in range(N_repeat):
        taux,Saux = Signal(Name,fs,rand_displace=1)        
        Sout = np.append(Sout,Saux)
        try:
            tout = np.append(tout,taux+tout[-1])
        except:
            tout = taux
    return tout,Sout





def Emitter_attitude(Name,time):
    if Name == 'Airplane':
        h = 10*km
        phi_o = (np.random.rand(1)-0.5)*np.pi*2 #random direction of first apearence
        Ro = 800*km # first distance
        if phi_o<np.pi/4 and phi_o>-np.pi/4: #coming from east
            direction = (np.random.rand(1)-0.5)*np.pi+np.pi #random direction    
        if phi_o<np.pi*3/4 and phi_o>np.pi/4: #coming from north
            direction = (np.random.rand(1)-0.5)*np.pi+np.pi*3/2 #random direction    
        if phi_o>np.pi*3/4 or phi_o<-np.pi*3/4: #coming from west
            direction = (np.random.rand(1)-0.5)*np.pi #random direction    
        if phi_o<-np.pi/4 and phi_o>-np.pi*3/4: # coming from south
            direction = (np.random.rand(1)-0.5)*np.pi+np.pi/2 #random direction    

        Xo = Ro*np.sin(phi_o)
        Yo = Ro*np.cos(phi_o)
        speed = 800*km_h
        X = Xo+np.sin(direction)*speed*time
        Y = Yo+np.cos(direction)*speed*time
        Z = h
    return X,Y,Z


def RFI_TX(Signals,powers,X,Y,Z):
    
#    Signals = matrix with the generated signal trains
#    powers = list with the Tx power for each signal type.
#    X,Y,Z = attitude of the emitter
     a=1
    
    
def RX_gain(X,Y,Z):
    #2 dimension sinc
    Go = 1e5
    R_half = 1.8954 # half power point of the sinc function
    beamwidth = 1*np.pi/180
    k = R_half*2/Go/np.tan(beamwidth/2)
    gain = np.sinc(np.sqrt(X**2+Y**2)/np.pi*k)*Go
    
    # returns the gain value in linea units
    return gain

def plot_rx_pattern():
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    from mpl_toolkits.mplot3d import Axes3D
    #2 dimension sinc
    lmax= 10000
    Go = 1e5
    R_half = 1.8954 # half power point of the sinc function
    beamwidth = 1*np.pi/180
    k = R_half*2/Go/np.tan(beamwidth/2)
    
    x = np.linspace(-lmax,lmax,100)
    y = np.linspace(-lmax,lmax,100)
    [X,Y] = np.meshgrid(x,y)
    Z = np.sinc(np.sqrt(X**2+Y**2)/np.pi*k)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    
    plt.show()

    
    
def Range(X,Y,h):
    R = np.sqrt(X**2+Y**2+h**2)
    return R

def Phase(R,fc):
    lda = 3e8/fc
    phase = 2*np.pi/lda*R
    return phase

def FreeSpaceLoss(R_m,fc_Hz):
    #receives Range in m and centre frequency in Hz    
    FSPL = 20*np.log10(R_m) + 20*np.log10(fc_Hz/MHz) - 27.55 
    return FSPL


#%% ADS-B
fs = 5000*MHz
[t1,S1] = Signal('ADS-B',fs)

plt.figure()
plt.plot(t1,S1)


[t1,S1] = gen_signal_train('ADS-B',fs,3)
plt.figure()
plt.plot(t1,S1)

N = len(S1)
P_fft = np.abs(np.fft.fftshift(np.fft.fft(S1)/len(S1)))**2/50
f_fft = np.fft.fftshift(np.fft.fftfreq(N, d=1/fs))

plt.figure()
plt.plot(f_fft,10*np.log10(P_fft*1e3))

# With noise
# Summing the noise to the signal train

sig = 1/10
mu = 0
Noise = Gauss_Noise(t1,fs,sig,mu)
S1 = S1 + Noise
plt.figure()
plt.plot(t1,S1)

    
    
#%% DME

[t1,S1] = Signal('DME',5000*MHz)
[t1,S2 ]= Signal('DME',5000*MHz)
[t1,S3 ]= Signal('DME',5000*MHz)

plt.figure()
plt.plot(t1,S1)
plt.plot(t1,S2)
plt.plot(t1,S3)


[t1,S1] = gen_signal_train('DME',5000*MHz,5)
plt.figure()
plt.plot(t1,S1)


#%% vuelos simulados avion.

t = np.linspace(0,3*hr,1000)
plt.figure()
k = np.linspace(0,1,1000)
X = 800*km*np.sin(2*np.pi*k)
Y = 800*km*np.cos(2*np.pi*k)
plt.plot(X,Y,'b')
for i in range(50):
    X,Y,Z = Emitter_attitude('Airplane',t)
    plt.plot(X,Y)
    plt.plot(X[0],Y[0],'ro')
    


