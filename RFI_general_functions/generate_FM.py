# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 09:11:26 2019

@author: f.divruno
"""
import numpy as np
import scipy.signal as signal

def generate_FM(audio,audio_SampleRate, Fc, SampleRate,freq_dev):
    N_audio = len(audio)
    t_audio = np.linspace(0,N_audio/audio_SampleRate,N_audio)
    N_FM = int(t_audio[-1]*SampleRate)
    t_FM = np.linspace(0,t_audio[-1],N_FM)
    audio_interp = np.interp(t_FM,t_audio,audio)
#    audio_resampled = signal.resample(audio,int(N_FM))/max(audio) #resampled and normalized to 1
    t_FM = np.linspace(0,t_audio[-1],N_FM)
    
    phase = np.cumsum(audio_interp)*2*np.pi*freq_dev/audio_SampleRate
    
    FM_signal = np.cos(2*np.pi*Fc*t_FM + phase)
    return FM_signal


#test
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    t_max = 1
    f_audio = 100
    audio_SampleRate = 44100 #Hz
    freq_dev = 2e3
    Fc = 1e6 
    Sample_rate = Fc*10
    
    t_audio = np.linspace(0,t_max,t_max*audio_SampleRate)
    audio = np.sin(2*np.pi*f_audio*t_audio)
    
    FM_signal = generate_FM(audio,audio_SampleRate,Fc,Sample_rate,freq_dev)
    
    
    plt.figure()
    f_fft = np.fft.fftfreq(len(FM_signal),d=1/Sample_rate)
    V_fft = np.abs(np.fft.fft(FM_signal))
    plt.plot(f_fft,20*np.log10(V_fft))
        