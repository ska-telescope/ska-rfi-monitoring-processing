# -*- coding: utf-8 -*-
"""
Created on Sun Jun 23 22:58:50 2019

@author: f.divruno
"""

#adcGain - compute adc gain
import numpy as np

OPTIMUM_THRES = .335 #sigma

def adcGain2(PdBm, nBits, adcVfs, thresOpt = OPTIMUM_THRES):
    """ PdBm power into ADC
        Nbit number of ADC bits
        adcVfs
    return attenuation required to optimitze input to
    optimum 4 bit correlator efficiency
    """
    Vrms = np.sqrt(10**(PdBm/10.)*.001*50.) #rms of input
    PdBmThresOpt = thresOpt * Vrms #Optimum threshold for the input
    adcThres = adcVfs/(2.**(nBits-1))
    Vratio = PdBmThresOpt/adcThres
    GdB = 20.*np.log10(Vratio)
    return(GdB)
    
if __name__ == '__main__':
    Pin = 10.
    nBits = 2.
    adcVfs = 1.0
    adcGain = adcGain2(Pin,nBits, adcVfs)
    print("To optimize ADC gain for 4 bit correlation add: ",  adcGain, " dBm of gain")
    
    calcThres = np.sqrt(10**((adcGain + Pin)/10.)*.001*50.)
    if np.abs((calcThres - OPTIMUM_THRES)/OPTIMUM_THRES > .001):
        print("adcGain Test Failed")
    else:
        print("adcGain Test Passed")

