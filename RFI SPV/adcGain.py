# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 16:15:33 2019

@author: f.divruno
"""

#adcGain - compute adc gain
import numpy as np

OPTIMUM_THRES = .335 #sigma

def adcGain(PdBm, nBits, adcVfs, thresOpt = OPTIMUM_THRES):
    """ PdBm power into ADC
        Nbit number of ADC bits
        adcVfs
    return attenuation required to optimitze input to
    optimum 4 bit correlator efficiency
    """
#    Vrms = np.sqrt(10**(PdBm/10.)*.001*50.) #rms of input
#    PdBmThresOpt = thresOpt * Vrms #Optimum threshold for the input
#    adcThres = adcVfs/(2.**(nBits-1))
#    Vratio = PdBmThresOpt/adcThres
#    GdB = 20.*np.log10(Vratio)

    Lev = adcVfs/2**(nBits-1)
    Vrms_i =  np.sqrt(10**(PdBm/10.)*.001*50.) #rms of input
    G = (Lev/Vrms_i/thresOpt)**2
    GdB = 10*np.log10(G)

    return(GdB)
    
if __name__ == '__main__':
    Pin = 10.
    nBits = 4.
    adcVfs = 1.0
    adcGain = adcGain(Pin,nBits, adcVfs)
    print("To optimize ADC gain for 4 bit correlation add: ",  adcGain, " dB of gain")
    
    Vadc_in = np.sqrt(10**((adcGain + Pin)/10.)*.001*50.) # RMS voltage
    Lev = adcVfs/2**(nBits-1)
    calcThres = Lev/Vadc_in
    
    if np.abs((calcThres - OPTIMUM_THRES)/OPTIMUM_THRES > .001):
        print("adcGain Test Failed")
    else:
        print("adcGain Test Passed")
