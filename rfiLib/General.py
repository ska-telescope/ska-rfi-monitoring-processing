# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 00:04:14 2019
    Collection of functions 

@author: f.divruno
"""

import scipy.io as sio
#import astropy.coordinates as Coord
import astropy.units as u
import numpy as np

def Coord_to_nparray(coord):
    
    A = np.array([coord[0].si.value,coord[1].si.value,coord[2].si.value])
    
    return A

def saveAntInData(antL, filePrefixName):
    """ Saves data to MATLAB file with vectors named 'time' & 'signal' """
    antDict = {}
    for ant in antL:
        antDict['height'] = ant.height
        antDict['Name'] = ant.Name
        antDict['lat'] = ant.lat.value
        antDict['lon'] = ant.lon.value
        antDict['Az'] = ant.Pointing['Azimuth'].value
        antDict['El'] = ant.Pointing['Elev'].value
        antDict['SampleRate'] = ant.SampleRate
        antDict['Duration'] = ant.Duration
        antDict['time'] = ant.time
        antDict['sky_source_rx'] = ant.sky_source_rx
        antDict['Rx_signal'] = ant.Rx_signal
        sio.savemat(filePrefixName + 'AntIn' + ant.Name.capitalize(), antDict)

def saveAdcInData(antL, filePrefixName):
    """ Saves data to MATLAB file with vectors named 'time' & 'signal' """
    antDict = {}
    for ant in antL:
        antDict['height'] = ant.height
        antDict['Name'] = ant.Name
        antDict['lat'] = ant.lat.value
        antDict['lon'] = ant.lon.value
        antDict['Az'] = ant.Pointing['Azimuth'].value
        antDict['El'] = ant.Pointing['Elev'].value
        antDict['SampleRate'] = ant.SampleRate
        antDict['Duration'] = ant.Duration
        antDict['time'] = ant.time
        antDict['ADC_input_rx'] = ant.ADC_input_rx
        antDict['ADC_input_sky'] = ant.ADC_input_sky

#        sio.savemat(filename+'_'+str(i), {"time":Telescope_list[i].time,"ADC_output_rx":Telescope_list[i].ADC_output_rx,"ADC_output_sky":Telescope_list[i].ADC_output_sky })
        # instead of saving the time vector save the sample rate, the time vector can be calculated.
        sio.savemat(filePrefixName + 'AdcIn' + ant.Name.capitalize(), antDict)

def saveAdcOutData(antL, filePrefixName):
    """ Saves data to MATLAB file with vectors named 'time' & 'signal' """
    antDict = {}
    for ant in antL:
        antDict['height'] = ant.height
        antDict['Name'] = ant.Name
        antDict['lat'] = ant.lat.value
        antDict['lon'] = ant.lon.value
        antDict['Az'] = ant.Pointing['Azimuth'].value
        antDict['El'] = ant.Pointing['Elev'].value
        antDict['SampleRate'] = ant.SampleRate
        antDict['Duration'] = ant.Duration
        antDict['time'] = ant.time
        antDict['ADC_input_rx'] = ant.ADC_output_rx
        antDict['ADC_input_sky'] = ant.ADC_output_sky

#        sio.savemat(filename+'_'+str(i), {"time":Telescope_list[i].time,"ADC_output_rx":Telescope_list[i].ADC_output_rx,"ADC_output_sky":Telescope_list[i].ADC_output_sky })
        # instead of saving the time vector save the sample rate, the time vector can be calculated.
        sio.savemat(filePrefixName + 'AdcOut' + ant.Name.capitalize(), antDict)
        
def loadAntInData(antL, filePrefixName):
    antDict = {}
    for ant in antL:
        antDict = sio.loadmat(filePrefixName + 'AntIn' + ant.Name.capitalize())
        ant.height = antDict['height'][0][0]
        ant.Name = antDict['Name'][0]
        ant.lat = antDict['lat'][0][0]*u.deg
        ant.lon = antDict['lon'][0][0]*u.deg
        ant.Pointing = dict(Elev=antDict['El'][0][0]*u.deg,Azimuth=antDict['Az'][0][0]*u.deg)
        ant.SampleRate = antDict['SampleRate'][0][0]
        ant.Duration = antDict['Duration'][0][0]
        ant.time = antDict['time'][0]
        ant.sky_source_rx = antDict['sky_source_rx'][0]
        ant.Rx_signal = antDict['Rx_signal'][0]
    return antL
        
