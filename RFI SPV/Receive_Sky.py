# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 23:52:03 2019

@author: f.divruno
"""

import General as General

'''----------------------------------

    ----------------------------------
'''

def Receive_Sky(Telescope_list, Sky_source_list, SamplingRate, Duration,plot_flag=0):
    
    max_delay_samples = (0) #Calculates the maximum delay in the signals to obtain same length
    total_samples = len(Telescope_list[0].Rx_signal)
    
    # will take a reference telescopewhere the received signal has zero delay.
    
    for i in range(len(Telescope_list)):
        print('\n\nTelescope: ' + Telescope_list[i].Name)
        
        for j in range(len(Sky_source_list)):
            print('\nSource: ' + Sky_source_list[j].Name)
            #for each emitter calculates the Range and the angle.
            Pos_Rx = General.Coord_to_nparray(Telescope_list[i].Pos)
            aux =  Sky_source_list[j].rx_sky(SamplingRate,Duration,Pos_Rx)
            Telescope_list[i].sky_source_rx = aux[0:total_samples] # here is only the sky source
            Telescope_list[i].Rx_signal += aux[0:total_samples]    # here is also the RFI 

    if plot_flag:
        Telescope_list[0].plot_signal('abs','Sky')
        Telescope_list[0].plot_signal('abslog','Sky')
        Telescope_list[0].plot_spectrum('abs','Sky')
    
    return Telescope_list


