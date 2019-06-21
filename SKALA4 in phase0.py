# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:28:07 2019

Reads a TDD data file from LFAA PHASE 0 experiment

@author: f.divruno
"""
import os, os.path
import numpy as np
import matplotlib.pyplot as plt
from math import * 
import matplotlib
import RFI_general_functions as RFI
from scipy import signal 
#from RFI_general_functions import * # RFI functions.

font = {'family' : 'DejaVu Sans','weight' : 'normal','size'   : 22}
matplotlib.rc('font', **font)


##%%
#
## A class that will downsample the data and recompute when zoomed.
#class DataDisplayDownsampler(object):
#    def __init__(self, xdata, ydata):
#        self.origYData = ydata
#        self.origXData = xdata
#        self.max_points = 1000
#        self.delta = xdata[-1] - xdata[0]
#
#    def downsample(self, xstart, xend):
#        # get the points in the view range
#        mask = (self.origXData >= xstart) & (self.origXData <= xend)
#        # dilate the mask by one to catch the points just outside
#        # of the view range to not truncate the line
#        mask = np.convolve([1, 1], mask, mode='same').astype(bool)
#        # sort out how many points to drop
#        ratio = max(np.sum(mask) // self.max_points, 1)
#
#        # mask data
#        xdata = self.origXData[mask]
#        ydata = self.origYData[mask]
#
#        # downsample xdata
#        xdata = xdata[::ratio]
#        # calculate max peak for y data
#        ydata = np.reshape(ydata,[len(ydata)//ratio,ratio])
#        ydata = np.max(ydata,1)
##        ydata = ydata[::ratio]
#
#        print("using {} of {} visible points".format(
#            len(ydata), np.sum(mask)))
#
#        return xdata, ydata
#
#    def update(self, ax):
#        # Update the line
#        lims = ax.viewLim
#        if np.abs(lims.width - self.delta) > 1e-8:
#            self.delta = lims.width
#            xstart, xend = lims.intervalx
#            self.line.set_data(*self.downsample(xstart, xend))
#            ax.figure.canvas.draw_idle()
#
#def plot_max_peak(xdata,ydata):
#    d = DataDisplayDownsampler(xdata, ydata)
#    fig, ax = plt.subplots()
#
#    # Hook up the line
##    d.line, = ax.plot(xdata, ydata)
#    d.line = ax.plot(xdata, ydata) 
#    ax.set_autoscale_on(False)  # Otherwise, infinite loop
#    
#    # Connect for changing the view limits
#    ax.callbacks.connect('xlim_changed', d.update)
#    ax.set_xlim(xdata.min(), xdata.max())
#    plt.show()
#    return d
#
#
#
#def read_phase0_data(directory,filenum):
#    files = os.listdir(directory)
#    # Each file is 1000 seconds, 4096 freq points, thats 16 minutes aprox
#    # there are aprox 4 files per hour.
#    # Data reduction:
#    # Each file look for max, min, 99%, 95%, 90%, limt line for RFI.
#    #
#        
#    header = 131072
#    if filenum == 'all':
#        N_files = np.size(files)
#    else:
#        N_files = filenum
#    time_vect = np.zeros([N_files])
#    data = np.zeros([N_files,header]).astype('int8')
#    i = 0
#    for i in range(N_files):
#        f = files[i]
#        fullpath = os.path.join(directory, f)
#        if os.path.splitext(fullpath)[1] == '.tdd':
#           with open(fullpath,'r') as f:
#               header = np.fromfile(f, dtype=np.dtype('>f8') , count = 1)#dtype=np.floating point 8 bytes
#               data[i,:] = np.fromfile(f, dtype=np.dtype('>i1') , count = 131072)#dtype=np.int64
#        sec = os.path.splitext(fullpath)[0][-2::]
#        mins = os.path.splitext(fullpath)[0][-4:-2]
#        hours = os.path.splitext(fullpath)[0][-6:-4]
#        time_vect[i] = int(hours)*60*60+int(mins)*60+int(sec)
#        i += 1
#        print(str(i))
#        
#    return time_vect,data
#
#
#def plot_waterfall(freq,data,title):
#      
#    fig = plt.figure(figsize=(20,12))
#    ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
#    #cbax = fig.add_axes((0.85, 0.05, 0.95, .95))
#    
#    
#    left = freq[0]
#    right = freq[-1]
#    bottom = 0
#    top = 2000
#    data_log = (10*np.log10(data))
#    cim = ax.imshow( data_log,
#        origin='lower', interpolation='nearest',
#        cmap= 'jet',
#        extent=(left, right, bottom, top),
#        )
#      
#    ax.set_aspect(abs(right-left) / abs(top-bottom))
#    plt.xlabel('MHz')
#    plt.ylabel('time')
#    plt.title(title)   
#
#def idle_fun(val):
#    return val
#
#def power_spectrum(data):
#    N = np.size(data,1) 
#    P = np.abs(data[:,0:int(N/2)].astype('float32'))**2
#    
#    return P
#
#def partition(signal,max_len, func1, func2= idle_fun):
#    N_files = np.size(signal,0)
#    N = np.size(signal,1)
##    data = np.zeros([0,N])
#    if N_files*N >  max_len: # steps to calculate all the FFTs
#        steps = int(N_files*N/max_len)
#        N_step = int(N_files/steps)
#        for i in range(steps):
#            A = func2(func1(signal[int(i*N_step):int((i+1)*N_step)]))
#            if i ==0:
#                data = A
#            else:
#                data = np.concatenate((data,A))
#            print(str(i+1)+' of '+str(steps))
#        if N_files > (i+1)*N_step:
#            A = func2(func1(signal[int((i+1)*N_step)::]))
#            data = np.concatenate((data,A))
#    else:
#        data= func2(func1(signal))
#    return data
#
#
#def fft_calc(signal, power_calc ,plot_figs):
#    # if power_calc = 1: Returns the power spectrum in real frequencies and the freq vector in MHz.
#    # if power_calc = 0: Returns the complex voltage spectrum and the freq vector in MHz.
#    
#    max_length =  8000*4096/10
#    N = np.size(signal,1)
#    N_files = np.size(signal,0)
#    
#    fs = 800e6
#    
#    freq = (np.fft.fftfreq(N, d=1/fs))
#    if power_calc ==1:
#        data_fft = partition(signal,max_length,np.fft.fft,power_spectrum)/N
#        freq = freq[0:int(N/2)]
#    else:
#        data_fft = partition(signal,max_length,np.fft.fft)/N
#    
#    plot_all = 0 # this is for debugging, if set to 1 generates 1 plot for each file.
#    if plot_all ==1:
##        N_files = 1 #for  debugging
#        plt.figure()
#        for i in range(N_files):
#            D = abs(data_fft)**2
#            plt.plot(freq/1e6,10*np.log10(D))    
#        plt.title('power spectrum' )
#
#    if plot_figs==1:
#        ave = np.average(abs(data_fft)**2,0)
#        plt.figure()
#        plt.plot(freq/1e6,10*np.log10(ave))    
#        plt.title('Average of all the captures')
#    return [freq/1e6,data_fft]
#
#
#def total_power(freq,data,fo,B,Tlow,Thigh,plot_flag):
##    fo = 160 #MHz
##    B = 1 #MHz
##    low_th = 0.4e1
##    high_th = 2e17
#
#    fmin = fo-B/2
#    fmax = fo+B/2
#    
#    fstep = freq[1] - freq[0]
#    
#    ind1 = int((fmin-freq[0])/fstep)
#    ind2 = int((fmax-freq[0])/fstep)
#    
#    total_power = np.sum(data[:,ind1:ind2]**2,1)
#    
#    mask = np.ones(len(total_power), dtype=bool)
#    for i in range(np.size(total_power)):
#        if (total_power[i] < Tlow or total_power[i] > Thigh):
#            mask[i] = False
#          
#    Data2 = total_power[mask]
#    
#    if plot_flag ==1:
#        plt.figure()
#        plt.plot(10*np.log10((Data2)))
#        plt.title('freq = ' + str(fo) + ' MHz, B = ' + str(B) + ' MHz')
#    return Data2

#% Read only one file to debug:
#
#directory = "C:\\Users\\f.divruno\\Dropbox (SKA)\\14- RFI environment\\01- Australia\\Phase-0\\2019-03-31\\DATA\\RX-02_SKALA-4.0\\Pol-X"
##directory = "G:\\Team Drives\\LowBridging Phase-0\\TPM\\DATA\\2019-03-31 (1)\\TRIGGER\\RX-02_SKALA-4.0\\Pol-X"
#    
#time,time_data = read_phase0_data(directory,filenum=1)
#[freq1,Ant1_fft] = fft_calc(time_data, 0,0)
#Power = Ant1_fft*np.conjugate(Ant1_fft)
#Autocorr = np.abs(np.fft.ifft(Power))
#
##[freq2,Ant2] = read_trigger_data(base_dir,date,antenna2,pol,1)
##[freq3,Ant3] = read_trigger_data(base_dir,date,antenna3,pol,1)
##[freq4,Ant4] = read_trigger_data(base_dir,date,antenna4,pol,1)
#
#Ant1_ave =  10*np.log10(np.average(abs(Power),0))
#
#
#plt.figure()
#plt.plot(freq1,Ant1_ave)
#plt.grid('on')
#plt.title('SPD1 plot')
#
#plot_all = 0
#if plot_all ==1:
#    
#    plt.figure()
#    plt.plot(np.transpose(Autocorr))
#    plt.grid('on')
#    plt.title('Auto Correlation')
#    
#    plt.figure()
#    plt.plot(np.transpose(time_data))
#    plt.grid('on')
#    plt.title('Time domain data')

#%% Read all the time files to calculate statistics in bands.
directory = "C:\\Users\\f.divruno\\Dropbox (SKA)\\14- RFI environment\\01- Australia\\Phase-0\\2019-03-31\\DATA\\RX-02_SKALA-4.0\\Pol-X"

read_data_from_files = 0
if  read_data_from_files ==1:
    time,td_data = RFI.read_phase0_data(directory,files_num='all')
    [freq,Ant1_pow] = RFI.fft_calc(td_data,fs=800e6, 0, 0)
    #del time_data
    np.savez_compressed(directory + '\\'+'Phase0_SKALA4_full_day', freq=freq, SKALA4_pow=Ant1_pow,time=time, td_data=td_data)
else:
    Aux = np.load('C:\\Users\\f.divruno\\Dropbox (SKA)\\14- RFI environment\\01- Australia\\Phase-0\\Phase0_SKALA4_full_day.npz')
    SKALA4_pow = Aux['SKALA4_pow'] #in V**2
    SKALA4_time = Aux['time'] # in seconds
    SKALA4_freq = Aux['freq'] # in MHz

load_td_data = 1
if load_td_data:
    Aux = np.load(r'C:\Users\F.Divruno\Dropbox (SKA)\14- RFI environment\01- Australia\Phase-0\2019-03-31\DATA\RX-02_SKALA-4.0\Pol-X\Phase0_full_day_raw.npz')
    td_data = Aux['td_data'] # time domain raw data
    

#%% filtering the td signal

fs = 800e6
fn = fs/2
f1 = 137.7e6/fn
f2 = 138e6/fn
ind = 49
b, a = signal.butter(3, [f1,f2],btype='bandpass',output='ba')

zi = signal.lfilter_zi(b, a)
z, _ = signal.lfilter(b, a, td_data[ind,:], zi=zi*td_data[ind,0])
#Apply the filter again, to have a result filtered at an order the same as filtfilt:

z2, _ = signal.lfilter(b, a, z, zi=zi*z[0])
#Use filtfilt to apply the filter:
y = signal.filtfilt(b, a, td_data[ind,:])

plt.figure()
plt.plot( td_data[ind,:], 'b', alpha=0.75)
plt.plot( z, 'r--', z2, 'r', y, 'k')
plt.legend(('noisy signal', 'lfilter, once', 'lfilter, twice','filtfilt'), loc='best')

[freq,P_orig] = RFI.fft_calc(td_data[ind,:],800e6, 1, 0)
[freq,P_filt] = RFI.fft_calc(y,800e6, 1, 0)

plt.figure()
plt.plot(freq,10*np.log10(P_orig),'r',freq,10*np.log10(P_filt),'b')

#%% maximum values

f= SKALA4_freq
data= SKALA4_pow
title= 'Max values - Phase0 - SKALA4.0 - 20190331'
RFI.plot_percentile(f,data,100,title='Maximum values - '+title)


#%% percentiles
perc = 100
f= SKALA4_freq
data= SKALA4_pow
title= 'Phase0 - SKALA4.0 - 20190331'
RFI.plot_percentile(f,data,perc,title)

#%% Spectrogram
Fstart = 0
Fstop = 30
time = SKALA4_time
f= SKALA4_freq
data= SKALA4_pow
title= 'Phase0 - SKALA4.0 - 20190331'

RFI.plot_spectrogram(time/3600,f,data,'Spectrogram -'+title,Fstart,Fstop,Tstart=0,Tstop=0)


#%% total power in a band
# This allows to see what is the time domain behaviour of the interference in this band.
time = SKALA4_time
f= SKALA4_freq
data= SKALA4_pow
title= 'Phase0 - SKALA4.0 - 20190331'

fo = 200 #MHz  # Aeronautic Comms band
B = 300 #MHz

#fo = 128 #MHz  # Aeronautic Comms band
#B = 16 #MHz

#fo = 138.5 #MHz #This is Orbcom
#B = 2 #MHz


#fo = 148.5 #MHz VHF radio?
#B = 1 #MHz

#fo = 255 #MHz Milstar satellite
#B = 30 #MHz

fo = 112.5 #MHz Milstar satellite
B = 9 #MHz

fo = 16 #MHz  # low freq stuf
B = 25 #MHz



low_th = 0
high_th = 1e200

time2,Data2 = RFI.total_power(time*3600,f,data,fo,B,low_th,high_th,0)

plt.figure()
plt.plot(time2/3600,10*np.log10(Data2))
plt.grid(True,'both') 
plt.title(title + ' - Fo=' + str(fo) + ' MHz - B=' + str(B) + 'MHz')

#%% Corss correlation of triggered samples
#base_dir = "G:\\Team Drives\\LowBridging Phase-0\\TPM\\DATA\\"
#date = "2019-03-15"
#
#ant1 = "RX-06_SKALA-4.1"
#ant2 = "RX-02_SKALA-4.0"
#ant3 = "RX-05_SKALA-2"
#ant4 = "RX-01_EDA-2"
#
#
#antenna1 = ant1
#antenna2 = ant2
#antenna3 = ant3
#pol1 = "Pol-X"
#pol2 = 'Pol-X'
#
#
#dir1 = base_dir + date + "\\TRIGGER\\" + antenna1 + "\\"+ pol1+"\\"
#dir2 = base_dir + date + "\\TRIGGER\\" + antenna2 + "\\"+ pol2+"\\"
#dir3 = base_dir + date + "\\TRIGGER\\" + antenna3 + "\\"+ pol1+"\\"
#
#
#filepos = 0
#time_data1 = read_phase0_data(dir1,filepos)
#time_data2 = read_phase0_data(dir2,filepos)
##time_data3 = read_phase0_data(dir3,filepos)
#[freq1,Ant1_fft] = fft_calc(time_data1, 0)
#[freq2,Ant2_fft] = fft_calc(time_data2, 0)
#Crosscorr = np.abs(np.fft.ifft(Ant1_fft*np.conjugate(Ant2_fft)))
#
#
##[freq2,Ant2] = read_trigger_data(base_dir,date,antenna2,pol,1)
##[freq3,Ant3] = read_trigger_data(base_dir,date,antenna3,pol,1)
##[freq4,Ant4] = read_trigger_data(base_dir,date,antenna4,pol,1)
#
#plt.figure()
#plt.plot(np.reshape(Crosscorr,[131072]))
#plt.grid('on')
#plt.title('Corss Correlation')
#
#plt.figure()
#plt.plot(np.reshape(time_data1,[131072]))
#plt.plot(np.reshape(time_data2,[131072]))
#plt.grid('on')
#plt.title('Time domain data')
#plt.legend([antenna1+pol1,antenna2+pol2])
#
#plt.figure()
#plt.plot(freq1/1e6,10*np.log10(np.abs(np.reshape(Ant1_fft,[131072])**2)))
#plt.plot(freq2/1e6,10*np.log10(np.abs(np.reshape(Ant2_fft,[131072])**2)))
#plt.grid('on')
#plt.title('Freq domain')
#plt.legend([antenna1+pol1,antenna2+pol2])
