# -*- coding: utf-8 -*-
"""
Created on Wed May 22 02:42:48 2019

@author: F.Divruno
"""
import matplotlib.pyplot as plt
import RFI_general_functions as RFI
import numpy as np
import sys
import os, os.path

ubuntu =0

print(sys.argv[0])
if len(sys.argv) < 1:
    print('Need parameter Ubuntu =1 , Windows=0')
    sys.exit(3)
else:
    ubuntu = sys.argv[1]


if ubuntu:
# in ubuntu
    outdir = r'/mnt/data/MRO_rfidata/output/'
    indir = r'/mnt/data/MRO_rfidata/'
else:
# in windows
    indir = 'C:\\Users\\F.Divruno\\Dropbox (SKA)\\14- RFI environment\\01- Australia\\rfidata_mro\\'
    outdir = 'C:\\Users\\F.Divruno\\Dropbox (SKA)\\14- RFI environment\\01- Australia\\rfidata_mro\\results\\'


def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype='float32')
    ret[n:] = ret[n:] - ret[:-n]
    ret[0:n-1] = a[0:n-1]*n #try
    #return ret[n - 1:] / n
    return ret/ n



read_files = 0
if read_files:
    print('reading files...')
    [freq,data] = RFI.read_MRO_data(indir,outdir)
    np.savez_compressed(outdir + r'MRO_rfidata_full', freq=freq, data=data)
    
else:
    try:
        Aux = np.load(outdir+r'MRO_rfidata_full.npz')
        data = Aux['data']/10-107 #in V**2 originally, scaled to get to dBm
        freq = Aux['freq'] # in MHz   
    except:
        files = os.listdir(outdir)
#        print(files)
        data = np.zeros([0,29801]).astype('float32')
        for i in range(len(files)): # comment for debug
#        for i in range(6):
            if os.path.splitext(files[i])[1] == '.npz':
                print(files[i])
                Aux = np.load(outdir+files[i])
                data = np.concatenate((data,Aux.get('data')/10-107),0) #in V**2 originally, scaled to get to dBm
                freq = Aux['freq'] # in MHz   
        np.savez_compressed(outdir + r'MRO_rfidata_full', freq=freq, data=data)

#%% USe the data from 0 to 500 MHz
fmin = 0
fmax = 500
D = data[:,(freq>=fmin) & (freq<=fmax)]

ave = np.average(D,0)            
plt.figure()
plt.plot(freq[(freq>=fmin) & (freq<=fmax)],ave)
plt.savefig(outdir+ 'Average_all', dpi=100, bbox_inches='tight')


#plot percentiles 
title = 'MRO data'
perc = 100
RFI.plot_percentile(freq[(freq>=fmin) & (freq<=fmax)],D,perc,'dBm',title)
title = title +'-'+ str(perc)+' percentile'
plt.savefig(outdir+title, dpi=100, bbox_inches='tight')

title = 'MRO data'
perc = 90
RFI.plot_percentile(freq[(freq>=fmin) & (freq<=fmax)],D,perc,'dBm',title)
title = title +'-'+ str(perc)+' percentile'
plt.savefig(outdir+title, dpi=100, bbox_inches='tight')



#%% polyfit filter.
#
#b = np.array(data)
#d = np.array(data)
#c = np.array(data)
##c = np.zeros([np.size(a,0),np.size(a,1)-(n_ave-1)])
#
#for i in range(len(data)):
#    blfit = 3
#    if blfit != 0:
#        # subtract an Nth order polynomial to remove baseline fluctuations
#        coeffs = np.polyfit(freq,data[i],blfit)
#        fit = np.poly1d(coeffs)(freq)
##        print ("max spect residual: %0.1f" % np.sqrt(np.median((a[i]-fit)**2)))
#        b[i] = data[i] - fit
#        print (i)
#
##
##    ave = moving_average(a[i,:],n_ave)
###    c[i,:] = a[i,19:] - ave
##    c[i,:] = a[i,:] - ave
#    
#    miin = np.min(data,0)
#    d[i] = data[i]- miin
#    

#%%
    
#        
#plt.figure()
#plt.plot(f,np.transpose(a[-1,:]))
#plt.plot(f,fit)
#plt.xlabel('Freq')
#plt.ylabel('power dBm')
#
##plt.figure()
##plt.plot(np.transpose(a[-1,:]))
##plt.plot(ave)
#
#
#        
#plt.figure()
#plt.plot(f,np.transpose(b[1,:]))
#plt.plot(f,np.transpose(d[1,:]))
##plt.plot(f_ave,np.transpose(c[1,:]))
##plt.plot(f,np.transpose(c[1,:]))
#
##%%
#ave = np.average(b,0)
#plt.figure()
##plt.plot(f_ave,ave)
#plt.plot(f,ave)

#%% 
#RFI.plot_percentile(f,d,100,'dBm','MRO data, normalized with the minimum [dB]')
#title = 'MRO data normalized with polyfit'
#perc = 100
#RFI.plot_percentile(freq,b,perc,'dBm',title)
#title = title +'-'+ str(perc)+' percentile'
#plt.savefig(outdir+ title, dpi=100, bbox_inches='tight')
#
#title = 'MRO data'
#perc = 100
#RFI.plot_percentile(freq,data,perc,'dBm',title)
#title = title +'-'+ str(perc)+' percentile'
#plt.savefig(outdir+title, dpi=100, bbox_inches='tight')
#
#title = 'MRO data'
#perc = 90
#RFI.plot_percentile(freq,data,perc,'dBm',title)
#title = title +'-'+ str(perc)+' percentile'
#plt.savefig(outdir+title, dpi=100, bbox_inches='tight')
#
#



