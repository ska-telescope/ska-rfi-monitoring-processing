# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 12:26:56 2019

@author: f.divruno
"""

import numpy as np
import matplotlib.pyplot as plt

def plot_percentile(freq,data,percentile,outdir,dataunits='dBm',title=''):
    #the data is m x n and the m dimension will be used to calculate the percentile, so 
    #this must be the time dimension.
#    data_perc = np.percentile(data,percentile,axis=0)
    data_sorted = np.sort(data)
    index = int((len(data)-1)*percentile/100)
    data_perc = data_sorted[index] 
    plt.figure()
    if dataunits == 'dBm':
        plt.plot(freq,(data_perc))
    else:
        plt.plot(freq,10*np.log10(data_perc))
    plt.grid(True,'both')
    plt.title(str(percentile) + ' percentile - ' + title)
    plt.savefig(outdir + str(percentile) + ' percentile - ' + title, dpi=600, bbox_inches='tight')
    return data_perc
    