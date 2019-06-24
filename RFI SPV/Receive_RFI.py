# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 23:51:32 2019

@author: f.divruno
"""
import numpy as np

from Path import Path

ms = 1e-3
us = 1e-6
MHz = 1e6
GHz = 1e9
km = 1e3
minute = 60
hr = 60*minute
km_h = km/hr
k_bolt = 1.23e-38


def Receive_RFI(Telescope_list, Emitter_list,Duration,SamplingRate,plot_flag=0):
    max_delay_samples = (0) #Calculates the maximum delay in the signals to obtain same length
    total_samples = int(Duration*SamplingRate)
    
    for i in range(len(Telescope_list)):
        print('\n\nTelescope: ' + Telescope_list[i].Name)
        Telescope_list[i].Rx_signal = np.zeros(total_samples)
        Telescope_list[i].time = np.linspace(0,Duration,total_samples)
        
        for j in range(len(Emitter_list)):
            print('\nEmitter: ' + Emitter_list[j].Name)
            #for each emitter calculates the Range and the angle.
            for k in range(len(Emitter_list[j].Signals)):
                print('Signal: ' + Emitter_list[j].Signals[k].Name)
                fc = Emitter_list[j].Signals[k].CenterFreq
                Pos_Tx = np.array(Emitter_list[j].Pos)
                Pos_Rx = [Telescope_list[i].Pos[0].value,Telescope_list[i].Pos[1].value,Telescope_list[i].Pos[2].value]
                
                # Here something that calculates the range and the time delay.
                #  TO-DO: Path should take into account the relative movement of both source and receiver.
                FSPL,delay = Path(Pos_Tx[:,0],Pos_Rx,fc) # Note: delay only depends on the range.

                delay_samples = int(delay*SamplingRate)
                if delay_samples > max_delay_samples:
                    max_delay_samples = delay_samples
                    if max_delay_samples/SamplingRate > Duration:
                        print('Delay is %0.2f us, greater han UC duration' % (max_delay_samples/SamplingRate/us))
                        exit(3)
                        
                FSPL_times = 10**(-FSPL/20)
                
                #Attenuate: Antenna gain is not included at the moment
                Sig_aux = FSPL_times*Emitter_list[j].Signals[k].data
                
                #Delay:
                Sig_aux = np.roll(Sig_aux,-delay_samples) # the delay makes the signal arrive later to the receiver
                
                Telescope_list[i].Rx_signal[0:len(Sig_aux)] += Sig_aux

            
#    print('chopping the received signals with max delay')
    for i in range(len(Telescope_list)):
#        print('\n\nTelescope: ' + Telescope_list[i].Name)
        total_samples = int(Duration*SamplingRate)
        Telescope_list[i].time = Telescope_list[i].time[0:(total_samples-max_delay_samples)]
        Telescope_list[i].Rx_signal = Telescope_list[i].Rx_signal[0:int(total_samples-max_delay_samples)]

    if plot_flag:
        Telescope_list[0].plot_signal('abs','RFI')
        Telescope_list[0].plot_signal('abslog','RFI')
        Telescope_list[0].plot_spectrum('abs','RFI')
    
    return Telescope_list
