# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 01:28:59 2019
test RFI Use cases framework
Derived from f. divruno FrameWork.py and modified to use support functions in
a support library. Make sure this library is in your PYTHONPATH

@author: G. Hovey
@author: f.divruno
@revised: G. Hovey;added code to reference antenna positions from csv file using pandas
@revised: G. Hovey;organized and documented test parameters 
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import pandas as pd

#RFI support functions
"""
import rfiLib.Sky as Sky
from rfiLib import Emitter as Emitter
import rfiLib.Receiver as Receiver
import rfiLib.Receive_RFI as Receive_RFI
import rfiLib.Receive_Sky as Receive_Sky
import rfiLib.Apply_DISH as Apply_DISH
"""

from rfiLib.Sky import Sky
from rfiLib.Emitter import Emitter
from rfiLib.Receiver import Receiver
from rfiLib.Receive_RFI import Receive_RFI
from rfiLib.Receive_Sky import Receive_Sky
from rfiLib.Apply_DISH import Apply_DISH
from rfiLib.General import saveAntInData, loadAntInData
from rfiLib.General import saveAdcInData, loadAdcInData
from rfiLib.General import saveAdcOutData, loadAdcOutData

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

#file parameters
testCaseName = 'test1'
skaMidAntPosFileSpec = './skaMidAntPositions.csv'
randomSeed = 25.

#antenna pair to test
tstAnt1Key = 'SKA001'
tstAnt2Key = 'SKA005'

#antenna pointing az and el
antAzEl = dict(Elev=90*u.deg,Azimuth=0*u.deg)

#Receiver and temporal parameters
Band = 'B2'
Duration = 2.*ms
SamplingRate = 4*GHz # THis is the analog sampling rate

#ADC scaling
scaling = 'Correlator_opimized'

#Test configuration parameters
promptFlg = False #interactive mode prompts user at various processing steps
runFlg = True #can be used to skip over processing
saveFlg = True #results are saved if true
loadFlg = False #results are loaded if true
plot_signal = True #plot time series signal
plot_spectrum = True #plot spectrum
plot_corr = False   #plot correlation
          

#%% Generation of the test case

def prompt(promptStr='Press enter to continue, or anything else to abort'):
    if promptFlg:
        return input(promptStr)
    else:
        return ''

skaMidAntPos = pd.read_csv(skaMidAntPosFileSpec, comment='#', index_col=0)

#Generate the RFI sources or emitters:
if((prompt('Generate RFI Sources [enter]?')=='') & runFlg):
    rfiSrcL = list([])
    rfiSrcL.append(Emitter('rfiSrc1','Airplane',dict(height_i = 10*u.km, lat_i = -30*u.deg, lon_i=20*u.deg), Duration, SamplingRate,[],random_seed=randomSeed,forceSignals=1))
    rfiSrcL.append(Emitter('rfiSrc2','Airplane',dict(height_i = 10*u.km, lat_i = -30.44*u.deg, lon_i=19.5*u.deg), Duration, SamplingRate,[],random_seed=randomSeed,forceSignals=1))


    print('Created RFI sources: ')
    for a in rfiSrcL: 
        print(a.Name + ' - ' +  a.Emit_type)
    
else:
    raise SystemExit


#Generate the antenna receivers:
if((prompt('Generate Antenna & Receivers [enter]?')=='') & runFlg):
    antRxL = list()
    antRxL.append(Receiver(skaMidAntPos.loc[tstAnt1Key].name,
                    dict(Latitude = skaMidAntPos.loc[tstAnt1Key].lat*u.deg, 
                        Longitude = skaMidAntPos.loc[tstAnt1Key].lon*u.deg),
                        antAzEl, Duration, SamplingRate))
                    
    antRxL.append(Receiver(skaMidAntPos.loc[tstAnt2Key].name,
                    dict(Latitude = skaMidAntPos.loc[tstAnt2Key].lat*u.deg, 
                         Longitude = skaMidAntPos.loc[tstAnt2Key].lon*u.deg),
                         antAzEl, Duration, SamplingRate))


    print('Created antennas: ')
    for a in antRxL: 
        print(a.Name)
                    
else:
    raise SystemExit
 
#Generate the sky signal sources
if((prompt('Generate Sky sources [enter]?')=='') & runFlg):
    skySrcL = list()
    
    skySrcL.append(Sky('Sky_source1', dict(lat= -31.340773 *u.deg,lon= 21.44*u.deg), SamplingRate, Duration, Temperature = 10))
#    skySrcL = list([skySrc1])
    print('Created sky sources: ')
    for a in skySrcL: 
        print(a.Name)
    
else:
    raise SystemExit


#Calculate the received RFI at the antenna aperture
if((prompt('Compute RFI at antenna aperture [enter]?')=='') & runFlg):
    antRxL = Receive_RFI(antRxL, rfiSrcL,Duration,SamplingRate,plot_flag=0)
    # Each of the receivers has the data in Receiver.Rx_signal 
    print('computing RFI sources at antenna aperture')
else:
    raise SystemExit


#Calculate the received sky signal at the antenna aperture and sum it with the RFI
#The received signal is stored in Rx_signal, the sky signal only is stored in sky_signal_rx
if((prompt('Compute Sky source at antenna aperture [enter]?')=='') & runFlg):
    antRxL = Receive_Sky(antRxL,skySrcL, SamplingRate, Duration,plot_flag=0)
    print('computing Sky sources at antenna aperture')
else:
    raise SystemExit



#Apply signals at aperture to the receiver chain and save points along the chain
#The output signal is stored in Receiver.ADC_output_rx (with RFI) or ADC_output_sky (without RFI)
if((prompt('Apply aperture signals to rx chain model [enter]?')=='') & runFlg):
    antRxL = Apply_DISH(antRxL,Band,scaling, atten = 0) 
    print('Taking antenna aperture singal applying to analog signal chain of :')
    for a in antRxL: 
        print(a.Name)
else:
    raise SystemExit


    # Save the data to files.
if saveFlg:
    saveAntInData(antRxL, testCaseName)
    saveAdcInData(antRxL, testCaseName)
    saveAdcOutData(antRxL, testCaseName)
    

if loadFlg:
    antRxL2 = list([Receiver('ant1'),Receiver('ant2')])
    antRxL2 = loadAntInData(antRxL2, 'test1')
    

#Plot the results

if plot_signal:
    for antRx in antRxL:
        antRx.plot_signal('antIn','RFI','volt') #mode= absVolt, volt, powerLin, powerLog
        antRx.plot_signal('adcIn','RFI','volt')
        antRx.plot_signal('adcOut','RFI','volt') 
        


if plot_spectrum:
    for antRx in antRxL:
        antRx.plot_spectrum('antIn','RFI','power')
        antRx.plot_spectrum('adcIn','RFI','power')
        antRx.plot_spectrum('adcOut','RFI','power')


#%% Verification of the results:
        
'''
TO-DO: 
   Calculate the SNR at the input and different taps of the SC model
   Calculate the SNR degradation (with and without RFI)

'''   

#%% Calculate Correlation

if plot_corr:
#    Corr = abs(np.fft.ifft(np.fft.fft(antRxL[0].Rx_signal)*np.conjugate(np.fft.fft(antRxL[1].Rx_signal))))
    Corr = abs(np.fft.ifft(np.fft.fft(antRxL[0].Rx_signal)*np.conjugate(np.fft.fft(antRxL[1].ADC_output_rx))))
    plt.figure()
    plt.plot(Corr)
    plt.title('Correlation of RFI + signal')
 
#    Corr = abs(np.fft.ifft(np.fft.fft(antRxL[0].sky_source_rx)*np.conjugate(np.fft.fft(antRxL[1].sky_source_rx))))
    Corr = abs(np.fft.ifft(np.fft.fft(antRxL[0].Rx_signal)*np.conjugate(np.fft.fft(antRxL[1].ADC_output_sky))))
    plt.figure()
    plt.plot(Corr)
    plt.title('Correlation of intended signal')
