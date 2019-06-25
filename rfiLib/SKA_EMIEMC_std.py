# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 09:03:40 2019

@author: f.divruno
"""
import numpy as np
import matplotlib.pyplot as plt


def SKA_EMIEMC_std(fstart,fstop,N,plot):
    """ Create the EMI/EMC standard lines
    inputs:
        Fstart: starting frequency in MHz
        Fstop: stop frequency in MHz
        N: Number of points to generate
        plot: 0 - no plot
              1 - plot PSD dBm/Hz
              2 - plot received power dBm
    """
        
    
    freq = np.logspace(np.log10(fstart),np.log10(fstop),N) #[MHz]
    # Continuum limit line in PSD
    SKA_Continuum_PSD =-17*np.log10(freq[freq<2000])-192 #[dBm/Hz]
    SKA_Continuum_PSD = np.append(SKA_Continuum_PSD,-249*np.ones(np.size(freq[freq>=2000])))
    RBW_Continuum = freq*1e6*1/100

    # Spectral line treshold
    SKA_Line_PSD = -17*np.log10(freq[freq<2000])-177 #[dBm/Hz]
    SKA_Line_PSD = np.append(SKA_Line_PSD,-234*np.ones(np.size(freq[freq>=2000])))
    RBW_Line = freq*1e6*0.001/100


    SKA_Continuum_dBm = SKA_Continuum_PSD + 10*np.log10(RBW_Continuum)
    SKA_Line_dBm = SKA_Line_PSD + 10*np.log10(RBW_Line)

    if plot == 1:
        plt.figure(figsize=(10, 7), dpi=80, facecolor='w', edgecolor='k')
        plt.plot(freq, SKA_Continuum_PSD, label = 'Continuum PSD')
        plt.xlabel("frequency [MHz]")
        plt.ylabel("[dBm/Hz]")
        plt.xscale('log')
        plt.plot(freq,SKA_Line_PSD, label = 'Spectral Line PSD')
        plt.grid(True,'both')
        plt.title('SKA PSD Limits [dBm/Hz]')
        plt.legend()

    if plot == 2:
        plt.figure(figsize=(10, 7), dpi=80, facecolor='w', edgecolor='k')
        plt.plot(freq, SKA_Continuum_dBm, label = 'Continuum Rx power')
        plt.xlabel("frequency [MHz]")
        plt.ylabel("[dBm]")
        plt.xscale('log')
        plt.plot(freq,SKA_Line_dBm, label = 'Spectral Line Rx power')
        plt.grid(True,'both')
        plt.title('SKA PSD Limits [dBm]')
        plt.legend()

        
#   return [freq,SKA_Continuum_PSD, RBW_Continuum, SKA_Continuum_dBm, SKA_Line_PSD,RBW_Line, SKA_Line_dBm]
    return {'freq' : freq*1e6,'SKA_Continuum_PSD':SKA_Continuum_PSD,'RBW_Continuum': RBW_Continuum, 'SKA_Continuum_dBm':SKA_Continuum_dBm, 'SKA_Line_PSD':SKA_Line_PSD,'RBW_Line':RBW_Line, 'SKA_Line_dBm':SKA_Line_dBm}