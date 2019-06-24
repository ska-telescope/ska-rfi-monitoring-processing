# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 12:54:59 2019

@author: f.divruno
"""
import numpy as np
import matplotlib.pyplot as plt
from adcGain import adcGain

def sample(signal, f_s, f_s_new):
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


def quantize(signal, nbits, Vfs):
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



nbits = 5
Vfs = 1

t = np.linspace(0,1,200)
signal = np.sin(2*np.pi*t)*1

#PdBm, nBits, adcVfs, thresOpt = OPTIMUM_THRES):
Vrms = np.max(signal)/np.sqrt(2)
P = Vrms**2/50
PdBm = 10*np.log10(P/1e-3)

gain_V = 10.**(adcGain(PdBm,nbits,Vfs)/20.)


plt.figure()
plt.plot(signal*gain_V)
plt.title('analog amplified signal')

qsignal = quantize(signal*gain_V, nbits , Vfs)

plt.figure(figsize=[15,10])
plt.plot(qsignal)
plt.title('quantized and amplified ')
