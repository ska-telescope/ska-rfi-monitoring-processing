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
RFI Test case #1:
    Two airplanes transmitting ADS-B and DME signals.
    Two SKA-MID telescopes
    
    ----------------------------------
'''


Duration = 2*ms
SamplingRate = 3*GHz # THis is the analog sampling rate

#GEnerate two signals: for debug, signals are included in the emitter.
#    S  = Signal('DME',Duration,SamplingRate,1100*MHz,30,0)
#    S1 = Signal('ADS-B',Duration,SamplingRate,1090*MHz,30,0)
#    plt.figure()
#    plt.plot(S.time/us,S.data)
#    plt.plot(S.time/us,S1.data)    


#Generate the emitters:
print('Generating Emitters')
Aeroplane1 = Emitter('Aeroplane1','Aeroplane',dict(height_i = 10*u.km, lat_i = -30*u.deg, lon_i=20*u.deg), Duration, SamplingRate,[])
Aeroplane2 = Emitter('Aeroplane2','Aeroplane',dict(height_i = 10*u.km, lat_i = -30.44*u.deg, lon_i=19.5*u.deg), Duration, SamplingRate,[])
  
Emitter_list = list([Aeroplane1,Aeroplane2])

#Generate the receivers:
print('Generating Receivers')
SKA_MID1 = Receiver('SKA001',dict(Latitude = -30.71329*u.deg, Longitude = 21.449412*u.deg),dict(Elev=90*u.deg,Azimuth=0*u.deg), Duration, SamplingRate)
SKA_MID2 = Receiver('SKA132',dict(Latitude = -31.340773*u.deg, Longitude = 21.267823*u.deg),dict(Elev=90*u.deg,Azimuth=0*u.deg), Duration, SamplingRate)

Telescope_list = list([SKA_MID1,SKA_MID2]) # List with receiver without receiving data

 
#Calculate the received RFI
Telescope_list = Receive_RFI(Telescope_list, Emitter_list,Duration,SamplingRate,plot_flag=0)
# Each of the receivers has the data in Receiver.Rx_signal 

#Generate the sky signal
Sky_source = Sky('Sky_source1', dict(lat= -31.340773 *u.deg,lon= 21.44*u.deg), SamplingRate, Duration, Temperature = 1000)

Sky_source_list = list([Sky_source])

#Calculate the received sky signal and sum it to the RFI
Telescope_list = Receive_Sky(Telescope_list,Sky_source_list, SamplingRate, Duration,plot_flag=0)


#The received signal is stored in Rx_signal, the sky signal only is stored in sky_signal_rx
Telescope_list = Apply_DISH(Telescope_list)






#Plot the result
plot_signal = 0
if plot_signal:
    for i in range(len(Telescope_list)):
        Telescope_list[i].plot_signal('abs','RFI')


   
#Calculate Correlation
plot_corr = 1
if plot_corr:
    Corr = abs(np.fft.ifft(np.fft.fft(Telescope_list[0].Rx_signal)*np.conjugate(np.fft.fft(Telescope_list[1].Rx_signal))))
    plt.figure()
    plt.plot(Corr)
    plt.title('Correlation of RFI + signal')
 
    Corr = abs(np.fft.ifft(np.fft.fft(Telescope_list[0].sky_source_rx)*np.conjugate(np.fft.fft(Telescope_list[1].sky_source_rx))))
    plt.figure()
    plt.plot(Corr)
    plt.title('Correlation of intended signal')

