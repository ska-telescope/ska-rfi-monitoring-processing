# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 23:51:32 2019

@author: f.divruno
@revised: G.Hovey, 13 July 2019; changed exit to raise exception with more meanful message
"""
import numpy as np

from rfiLib.Path import Path

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
    total_samples = int(Duration*SamplingRate)
    delay_samples = np.zeros([len(Telescope_list),len(Emitter_list)])
    FSPL = np.zeros([len(Telescope_list),len(Emitter_list)])
    
    for i in range(len(Telescope_list)):
        print('\n\nTelescope: ' + Telescope_list[i].Name)
        Telescope_list[i].Rx_signal = np.zeros(total_samples)
        Telescope_list[i].time = np.linspace(1/SamplingRate,Duration,total_samples)
        
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
                FSPL[i,j],delay = Path(Pos_Rx,Pos_Tx[:,0],fc) # Note: delay calculation updated
                delay_samples[i,j] = int(delay*SamplingRate)

                        
                FSPL_times = 10**(-FSPL/20)
                

    # apply the delay to each emitter on each receiver
    for j in range(len(Emitter_list)):
        delay_emitter = delay_samples[:,j] - delay_samples[:,j].min() # calculate the differential delays
        # if there is a receiver with 0 delay, it must be delayed with max_delay, the receiver with the largest delay
        # gets the zero delay to be able to chop different parts of the same emitter
        delay_emitter -= delay_emitter.max()
        for i in range(len(Telescope_list)):
            for k in range(len(Emitter_list[j].Signals)):
                #Attenuate: Antenna gain is not included at the moment
                Sig_aux = FSPL_times[i,j]*Emitter_list[j].Signals[k].data
                
                Sig_aux = np.roll(Sig_aux,int(delay_emitter[i])) # the delay makes the signal arrive earlier to the receiver                
                Telescope_list[i].Rx_signal += Sig_aux[0:total_samples]    
                
                #added to calculate the maximum level at the antenna input:
                print('Ant: %s - Emitter: %s - Max level at antenna input: %f uV'%(Telescope_list[i].Name,Emitter_list[j].Name,np.max(Sig_aux)*1e6))
                print('Ant: %s - Emitter: %s - Delay: %f us'%(Telescope_list[i].Name,Emitter_list[j].Name,np.abs(delay_emitter[i]/SamplingRate)*1e6))
            
##    print('chopping the received signals with max delay')
#    for i in range(len(Telescope_list)):
##        print('\n\nTelescope: ' + Telescope_list[i].Name)
#        total_samples = int(Duration*SamplingRate)
#        Telescope_list[i].time = Telescope_list[i].time[0:(total_samples-max_delay_samples)]
#        Telescope_list[i].Rx_signal = Telescope_list[i].Rx_signal[0:int(total_samples-max_delay_samples)]

    if plot_flag:
        Telescope_list[0].plot_signal('abs','RFI')
        Telescope_list[0].plot_signal('abslog','RFI')
        Telescope_list[0].plot_spectrum('abs','RFI')
    
    return Telescope_list


