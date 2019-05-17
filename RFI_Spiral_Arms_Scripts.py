''' Author      : Dr. A. J. Otto
    Company     : MESA Solutions (Pty) Ltd.
    Date        : 03 May 2016

    Description : Scripts to present RFI data from February 2016 campaign using
                  the MESA Product Solutions (Pty) Ltd RTA-3 (Real Time
                  Analyser) and the passive Alaris Omni A0190-02 antenna. 


                modified 17/05/2019 by FDV to be used with Python 3.6
                
                '''
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import signal, integrate

mpl.rcParams['font.size'] = 16


def get_calibrated_data(f):
    data = f['calibrated_spectrum']
    freq = np.array(f['freqs'])

    data_max = np.max(data[0:-1, :], axis=0)
    data_median = np.median(data[0:-1, :], axis=0)
    data90 = np.percentile(data[0:-1, :], 90., axis=0)
    data99 = np.percentile(data[0:-1, :], 99., axis=0)

    return freq, data, data_max, data90, data99, data_median


def plot_calibrated_data(freqs1, data_max1, data991, data901, data_median1,
                         low1, high1):
    plt.figure(figsize=(15, 10))
    plt.semilogx(freqs1[low1:high1]/1e6, data_max1[low1:high1], '-r')
    plt.semilogx(freqs1[low1:high1]/1e6, data991[low1:high1], '-g')
    plt.semilogx(freqs1[low1:high1]/1e6, data901[low1:high1], '-b')
    plt.semilogx(freqs1[low1:high1]/1e6, data_median1[low1:high1], '--m')
    plt.legend(['Maximum', '99 Percentile', '90 Percentile', '50 Percentile'],
               loc=0, fontsize=12)
    plt.grid(True, 'major')
    plt.grid(True, 'minor')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Power [dBm]')
    plt.title('Measured Background RFI')
    plt.show()


def freq_to_chan(freqs, f):
    return np.array(abs(freqs - f)).argmin()


def get_time_occupancy(freqs1, data1, chan_low1, chan_high1, title):
    df = freqs1[1] - freqs1[0]
    # 2D Median Filter
    D_t = 1
    D_f = 10e6
    filteredData = signal.medfilt2d(data1[:, chan_low1:chan_high1],
                                    [2*int(D_t/2) + 1, 2*int(D_f/df/2) + 1])

    ''' This will sum the number of times the measured data is 1 dB above
        the noise floor calculated with the 2D median filter '''
    tmp = sum(i > j for (i, j) in zip(data1[:, chan_low1:chan_high1],
                                      filteredData + 1.))

    to = []
    for i in tmp:
        to.append(float(i)/float(data1.shape[0])*100.)

    plt.figure(figsize=(15, 10))
    plt.plot(freqs1[chan_low1:chan_high1]/1e6, to)
    plt.grid(True, 'major')
    plt.axis([freqs1[chan_low1]/1e6, freqs1[chan_high1]/1e6, -5, 105])
    plt.ylabel('Time Occupancy [%]')
    plt.xlabel('Frequency [MHz]')
    plt.title(title)
    plt.show()


def get_waterfall_section(freqs, data,  data_max, low, high, title):
    chan_low = np.array(abs(freqs - low)).argmin()
    chan_high = np.array(abs(freqs - high)).argmin()

    data = data[:, chan_low:chan_high]

    data_min = np.min(data[0:-1, :], axis=0)
    data_max = data_max[chan_low:chan_high]

    x = freqs[chan_low:chan_high]/1e6
    y = np.arange(np.shape(data)[0])
    X, Y = np.meshgrid(x, y)
    Z = data

    plt.figure(figsize=(15, 10))
    vmin = round(min(data_min), 0) - 1.
    vmax = round(max(data_max), 0)
    cf = plt.contourf(X, Y, Z, np.arange(vmin, vmax, 0.1))
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Time [s]')
    plt.title(title)
    cbar = plt.colorbar(cf, ticks=np.arange(vmin, vmax, 2.))
    cbar.ax.set_ylabel('Power [dBm]')
    plt.tight_layout()
    plt.show()


def get_int_power(f, freqs, data, low, high):
    chan_low = np.array(abs(freqs - low)).argmin()
    chan_high = np.array(abs(freqs - high)).argmin()
    df = freqs[1] - freqs[0]
    freqs = freqs[chan_low:chan_high]
    data = data[0:-1, chan_low:chan_high]
    # 2D Median Filter
    D_t = 1
    D_f = 10e6
    print ("Median Filter Window Size: ")
    print ([2*int(D_t/2) + 1, 2*int(D_f/df/2) + 1])
    filteredData = signal.medfilt2d(data,
                                    [2*int(D_t/2) + 1, 2*int(D_f/df/2) + 1])
    print ("Filtered Data/Noise Shape [C]: ", filteredData.shape)
    print ("Filtered Data/Noise [C]:")
    print (filteredData)

    # Aggregate (mean, max, 99 percentile etc.)
    ''' MEAN '''
    newData, filteredData =\
        np.nanmean(data, axis=0), np.nanmean(filteredData, axis=0)
    subtitle = "Mean"

    ''' MAX '''
    # newData, filteredData =\
    #     np.nanmax(newData, axis=0), np.nanmax(filteredData, axis=0)
    # subtitle = "Max"

    print ("Aggregate (" + subtitle + "): ")
    print ("Aggregate Data Shape [S]: ", newData.shape)
    print ("Aggregate Filtered Data/Noise Shape [C]: ", filteredData.shape)

    # Convert from Logarithmic to Linear values:
    newdata, filtereddata = 10.**(newData/10.), 10.**(filteredData/10.)

    # Calculate residual RFI spectrum
    # (RFI = Measured Data - Noise (Filtered Data))
    resdata = newdata - filtereddata
    resdata[resdata < 0] = 0
    filtereddata = newdata - resdata

    newData, filteredData, resData =\
        10.*np.log10(np.abs(newdata)), 10.*np.log10(np.abs(filtereddata)),\
        10.*np.log10(np.abs(resdata))

    # Logarithmic Plots
    plt.figure(figsize=(15, 10))
    plt.suptitle('Background RFI Measurements')
    plt.subplot(2, 1, 1)
    plt.plot(freqs/1e6, newData, '-b')
    plt.plot(freqs/1e6, filteredData, 'or')
    plt.grid(True, 'major')
    plt.grid(True, 'minor')
    plt.ylabel('Measured Spectrum [dBm]')
    plt.subplot(2, 1, 2)
    plt.plot(freqs/1e6, resData, 'bo')
    plt.grid(True, 'major')
    plt.grid(True, 'minor')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Residual RFI Spectrum [dBm]')
    plt.show()
    return df, newdata, resdata, filtereddata

#%% 
''' Site Location:  '''
''' *************** '''
''' CORE SITE (M48) '''
''' *************** '''

''' BAND 0 '''
low = 350e6
high = 2000e6

''' This is the path where the .h5 data file is located '''
path = 'calibrated_data/Day1_Core/M48/'
f1 = h5py.File(path + 'M48_BackgroundRFI_0_1455698254_Calibrated.spec.h5',
               'r')

''' *************** '''
'''    READ DATA    '''
''' *************** '''
''' This function will read the calibrated amplitude and frequency data from
    the calibrated .h5 file. It will also calculate the max, 50th, 90th and
    99th percentiles of the measured data '''
freqs1, data1, data_max1, data901, data991, data_median1 =\
    get_calibrated_data(f1)

''' This function returns the index number where a frequeny [Hz] is found '''
chan_low1 = freq_to_chan(freqs1, low)
chan_high1 = freq_to_chan(freqs1, high)

''' *************** '''
'''    PLOT DATA    '''
''' *************** '''
''' This function will plot the data obtained from get_calibrated_data() '''
plot_calibrated_data(freqs1, data_max1, data901, data991, data_median1,
                     chan_low1, chan_high1)
''' *************** '''

#%%
''' *************** '''
''' TIME OCCUPANCY  '''
''' *************** '''
''' This function will calculate the Time Occupancy '''
title = "UHF Satcom"
low = 245e6
high = 270e6
chan_low1 = freq_to_chan(freqs1, low)
chan_high1 = freq_to_chan(freqs1, high)
get_time_occupancy(freqs1, data1, chan_low1, chan_high1, title)
''' *************** '''

#%%
''' *************** '''
''' WATERFALL PLOTS '''
''' *************** '''
''' This function will calculate the Waterfall Plot for the given section: '''
title = "Carnarvon SABC"
low1 = 622.5e6
high1 = 625.5e6
get_waterfall_section(freqs1, data1, data_max1, low1, high1, title)

#%%
''' **************** '''
''' INTEGRATED POWER '''
''' **************** '''
''' This function will calculate the total integrated power over a band '''
# Read data and get integrated RFI power over band of interest
low = 350e6
high = 750e6
df, newdata, resdata, filtereddata =\
    get_int_power(f1, freqs1, data1, low, high)

# Integrate RFI power over specific band
nd, rfi, n = [], [], []
band_freqs = [(0.9*low, 1.1*high)]
for band in band_freqs:
    f_lo = freqs1.searchsorted(band[0])
    f_hi = freqs1.searchsorted(band[1])
    # Integrate Measured Data
    nd.append(integrate.trapz(newdata[f_lo:f_hi]/df, dx=df, axis=0))
    # Integrate RFI
    rfi.append(integrate.trapz(resdata[f_lo:f_hi]/df, dx=df, axis=0))
    # Integrate Noise
    n.append(integrate.trapz(filtereddata[f_lo:f_hi]/df, dx=df, axis=0))

ND, RFI, N =\
    10.*np.log10(np.abs(nd)), 10*np.log10(np.abs(rfi)), 10*np.log10(np.abs(n))

# print ("Integrated Measured Spectrum Power: ", ND, " dBm"
print ("Integrated RFI Spectrum Power: ", RFI, " dBm ", rfi, " mW")
print ("Integrated Noise Power: ", N, " dBm ", n, " mW\n")
''' *************** '''
