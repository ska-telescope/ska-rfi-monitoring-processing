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
import astropy.constants as const
import pandas as pd

#RFI support functions

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
testCaseName = 'test2'
skaMidAntPosFileSpec = './skaMidAntPositions.csv'
randomSeed = int(22)
maxDelay = 1e-3 *u.s#

#antenna pair to test
tstAnt1Key = 'SKA008'
tstAnt2Key = 'SKA133'

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
saveFlg = False #results are saved if true
loadFlg = False #results are loaded if true
plot_signal = True #plot time series signal
plot_spectrum = True #plot spectrum
plot_corr = True   #plot correlation
          

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
    rfiSrcL.append(Emitter('rfiSrc1','Airplane',dict(height_i = 10*u.km, lat_i = -30*u.deg, lon_i=20*u.deg), Duration+maxDelay.value, SamplingRate,[],random_seed=randomSeed,forceSignals=0))
    rfiSrcL.append(Emitter('rfiSrc2','Airplane',dict(height_i = 10*u.km, lat_i = -30.44*u.deg, lon_i=19.5*u.deg), Duration, SamplingRate,[],random_seed=randomSeed*2,forceSignals=0))
    rfiSrcL.append(Emitter('rfiSrc3','Airplane',dict(height_i = 10*u.km, lat_i = -31.7*u.deg, lon_i=20.5*u.deg), Duration, SamplingRate,[],random_seed=randomSeed*3,forceSignals=1))
    rfiSrcL.append(Emitter('rfiSrc4','Airplane',dict(height_i = 10*u.km, lat_i = -31*u.deg, lon_i=21*u.deg), Duration, SamplingRate,[],random_seed=randomSeed*4,forceSignals=0))
    rfiSrcL.append(Emitter('rfiSrc5','Airplane',dict(height_i = 10*u.km, lat_i = -30.7*u.deg, lon_i=21.3*u.deg), Duration, SamplingRate,[],random_seed=randomSeed*5,forceSignals=0))
    rfiSrcL.append(Emitter('rfiSrc6','Airplane',dict(height_i = 10*u.km, lat_i = -30.4*u.deg, lon_i=21.5*u.deg), Duration, SamplingRate,[],random_seed=randomSeed*6,forceSignals=1))
    rfiSrcL.append(Emitter('rfiSrc7','Airplane',dict(height_i = 10*u.km, lat_i = -30*u.deg, lon_i=21.8*u.deg), Duration, SamplingRate,[],random_seed=randomSeed*7,forceSignals=0))


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
    
    skySrcL.append(Sky('Sky_source1', dict(elev=70 *u.deg,az= 0*u.deg),antRxL[0].Pos, SamplingRate, Duration+maxDelay.value, Temperature = 20, random_seed=int(randomSeed*10)))
#    skySrcL.append(Sky('Sky_source2', dict(elev=60 *u.deg,az= 0*u.deg),antRxL[0].Pos, SamplingRate, Duration+maxDelay.value, Temperature = 10, random_seed=int(randomSeed*20) ))
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
#        antRx.plot_signal('adcOut','sky','volt') 


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
   Rudimental flagging alghorithm in time domain
   
'''   

#%% Calculate Correlation

if plot_corr:
#    Corr = abs(np.fft.ifft(np.fft.fft(antRxL[0].Rx_signal)*np.conjugate(np.fft.fft(antRxL[1].Rx_signal))))
    Corr = abs(np.fft.ifft(np.fft.fft(antRxL[0].ADC_output_rx)*np.conjugate(np.fft.fft(antRxL[1].ADC_output_rx))))
    plt.figure()
    plt.plot(Corr)
    plt.title('Correlation of RFI + signal, antennas: %s - %s'%(antRxL[0].Name,antRxL[1].Name))
  
#    Corr = abs(np.fft.ifft(np.fft.fft(antRxL[0].sky_source_rx)*np.conjugate(np.fft.fft(antRxL[1].sky_source_rx))))
    Corr = abs(np.fft.ifft(np.fft.fft(antRxL[0].ADC_output_sky)*np.conjugate(np.fft.fft(antRxL[1].ADC_output_sky))))
    plt.figure()
    plt.plot(Corr)
    plt.title('Correlation of intended signal output ADC, antennas: %s - %s'%(antRxL[0].Name,antRxL[1].Name))

#    Corr = abs(np.fft.ifft(np.fft.fft(antRxL[0].sky_source_rx)*np.conjugate(np.fft.fft(antRxL[1].sky_source_rx))))
#    plt.figure()
#    plt.plot(Corr)
#    plt.title('Correlation of intended signal input ADC')
#



