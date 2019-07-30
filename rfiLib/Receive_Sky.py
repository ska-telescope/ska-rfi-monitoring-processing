# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 23:52:03 2019

@author: f.divruno
"""

import rfiLib.General as General
import numpy as np


'''----------------------------------

    ----------------------------------
'''

def Receive_Sky(Telescope_list, Sky_source_list, SamplingRate, Duration,plot_flag=0):
    
    total_samples = int(SamplingRate*Duration)
    
    # will take a reference telescopewhere the received signal has zero delay.
    delaySamples = np.zeros([len(Telescope_list),len(Sky_source_list)])
    for i in range(len(Telescope_list)):
        print('\n\nTelescope: ' + Telescope_list[i].Name)
        
        for j in range(len(Sky_source_list)):
            print('\nSource: ' + Sky_source_list[j].Name)
            #for each emitter calculates the Range and the angle.
            Pos_Rx = General.Coord_to_nparray(Telescope_list[i].Pos)
            delaySamples[i,j] =  Sky_source_list[j].rx_sky(SamplingRate,Duration,Pos_Rx)

   
   # apply the delay to each emitter on each receiver
    for j in range(len(Sky_source_list)):
        delayEmitter = delaySamples[:,j]-delaySamples[:,j].min() # calculate the differential delays
        # if there is a receiver with 0 delay, it must be delayed with max_delay, the receiver with the largest delay
        # gets the zero delay to be able to chop different parts of the same emitter
        delayEmitter -= delayEmitter.max()
        for i in range(len(Telescope_list)):
            #Attenuate: Antenna gain is not included at the moment
            Sig_aux = Sky_source_list[j].data

            Sig_aux = np.roll(Sig_aux,int(delayEmitter[i])) # the delay makes the signal arrive earlier to the receiver                
            Telescope_list[i].Rx_signal += Sig_aux[0:total_samples]    
            Telescope_list[i].sky_source_rx = Sig_aux[0:total_samples] # here is only the sky source        
            

    if plot_flag:
        Telescope_list[0].plot_signal('abs','Sky')
        Telescope_list[0].plot_signal('abslog','Sky')
        Telescope_list[0].plot_spectrum('abs','Sky')
    
    return Telescope_list


