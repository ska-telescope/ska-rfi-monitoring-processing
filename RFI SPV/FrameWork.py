# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 01:28:59 2019
RFI Use cases framework


@author: f.divruno
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from Sky import Sky
from Emitter import Emitter
from Receiver import Receiver
from Receive_RFI import Receive_RFI
from Receive_Sky import Receive_Sky
from Apply_DISH import Apply_DISH


ms = 1e-3
us = 1e-6
MHz = 1e6
GHz = 1e9
km = 1e3
minute = 60
hr = 60*minute
km_h = km/hr
k_bolt = 1.23e-38

   
    

        
'''----------------------------------
RFI Signal Performance Verification

RFI Test case #1:
    Two airplanes transmitting ADS-B and DME signals.
    Two SKA-MID telescopes
    Band 2
    ----------------------------------
'''

#%% Generation of the test case

Duration = 2*ms
SamplingRate = 4*GHz # THis is the analog sampling rate
Band = 'B2'
scaling = 'Correlator_opimized'



#Generate the emitters:
if(input('Generate Emitters?')==''):
    Aeroplane1 = Emitter('Aeroplane1','Aeroplane',dict(height_i = 10*u.km, lat_i = -30*u.deg, lon_i=20*u.deg), Duration, SamplingRate,[])
    Aeroplane2 = Emitter('Aeroplane2','Aeroplane',dict(height_i = 10*u.km, lat_i = -30.44*u.deg, lon_i=19.5*u.deg), Duration, SamplingRate,[])
      
    Emitter_list = list([Aeroplane1,Aeroplane2])
else:
    exit(3)



#Generate the receivers:
if(input('Generate Receivers?')==''):
    SKA_MID1 = Receiver('SKA001',dict(Latitude = -30.71329*u.deg, Longitude = 21.449412*u.deg),dict(Elev=90*u.deg,Azimuth=0*u.deg), Duration, SamplingRate)
    SKA_MID2 = Receiver('SKA132',dict(Latitude = -31.340773*u.deg, Longitude = 21.267823*u.deg),dict(Elev=90*u.deg,Azimuth=0*u.deg), Duration, SamplingRate)
    Telescope_list = list([SKA_MID1,SKA_MID2]) # List with receiver without receiving data
else:
    exit(3)
 
#Calculate the received RFI
if(input('Receive RFI?')==''):
    Telescope_list = Receive_RFI(Telescope_list, Emitter_list,Duration,SamplingRate,plot_flag=0)
    # Each of the receivers has the data in Receiver.Rx_signal 
else:
    exit(3)


#Generate the sky signal
if(input('Generate Sky signal?')==''):
    Sky_source = Sky('Sky_source1', dict(lat= -31.340773 *u.deg,lon= 21.44*u.deg), SamplingRate, Duration, Temperature = 10)

    Sky_source_list = list([Sky_source])
else:
    exit(3)

if(input('Receive Sky signal?')==''):
    #Calculate the received sky signal and sum it to the RFI
    Telescope_list = Receive_Sky(Telescope_list,Sky_source_list, SamplingRate, Duration,plot_flag=0)
    #The received signal is stored in Rx_signal, the sky signal only is stored in sky_signal_rx
else:
    exit(3)

if(input('Apply DISH simplified model?')==''):
    Telescope_list = Apply_DISH(Telescope_list,Band,scaling, atten = 0) 
    # The signal inputing to the ADC is in the variable Receiver.ADC_input_rx or ADC_input_sky
    # The output signal is stored in Receiver.ADC_output_rx (with RFI) or ADC_output_sky (without RFI)
else:
    exit(3)

#to do:
    # Save the ADC output to files.
    

#Plot the results
plot_signal = 1
if plot_signal:
    for i in range(len(Telescope_list)):
        Telescope_list[i].plot_signal('Ant_in','RFI','abs')
        Telescope_list[i].plot_signal('ADC_in','RFI','abs')
        Telescope_list[i].plot_signal('ADC_out','RFI','abs')
        

plot_spectrum = 1
if plot_signal:
    for i in range(len(Telescope_list)):
        Telescope_list[i].plot_spectrum('Ant_in','RFI','abs')
        Telescope_list[i].plot_spectrum('ADC_in','RFI','abs')
        Telescope_list[i].plot_spectrum('ADC_out','RFI','abs')


#%% Verification of the results:
        
'''
TO-DO: 
   Calculate the SNR at the input and different taps of the SC model
   Calculate the SNR degradation (with and without RFI)

'''   
#%% Calculate Correlation
plot_corr = 1
if plot_corr:
#    Corr = abs(np.fft.ifft(np.fft.fft(Telescope_list[0].Rx_signal)*np.conjugate(np.fft.fft(Telescope_list[1].Rx_signal))))
    Corr = abs(np.fft.ifft(np.fft.fft(Telescope_list[0].Rx_signal)*np.conjugate(np.fft.fft(Telescope_list[1].ADC_output_rx))))
    plt.figure()
    plt.plot(Corr)
    plt.title('Correlation of RFI + signal')
 
#    Corr = abs(np.fft.ifft(np.fft.fft(Telescope_list[0].sky_source_rx)*np.conjugate(np.fft.fft(Telescope_list[1].sky_source_rx))))
    Corr = abs(np.fft.ifft(np.fft.fft(Telescope_list[0].Rx_signal)*np.conjugate(np.fft.fft(Telescope_list[1].ADC_output_sky))))
    plt.figure()
    plt.plot(Corr)
    plt.title('Correlation of intended signal')

