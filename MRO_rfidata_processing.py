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

ubuntu = 1
#
#print(sys.argv[0])
#if len(sys.argv) <= 1:
#    print('Need parameter Ubuntu =1 , Windows=0')
#    sys.exit(3)
#else:
#    ubuntu = sys.argv[1]


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
        print('Loading the complete set of saved data...')
        Aux = np.load(outdir+r'MRO_rfidata_full.npz')
        data = Aux['data']/10-107 #in V**2 originally, scaled to get to dBm
        freq = Aux['freq'] # in MHz   
    except:
        print('Failed...')
        print('Loading the saved data by files one by one...')
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

#%%
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


#%% Occupancy in Balt's method:

print ("Calculating occupancy...")
occup = np.zeros(len(D[0]))

freqs = np.array(freq[(freq>=fmin) & (freq<=fmax)])
amps = np.array(D)
blfit = 3

for i in range(len(amps)):
  if blfit != 0:
    # subtract an Nth order polynomial to remove baseline fluctuations
    coeffs = np.polyfit(freqs,amps[i],blfit)
    fit = np.poly1d(coeffs)(freqs)
    #print ("max spect residual: %0.1f" % np.sqrt(np.median((mspect-fit)**2))
    amps[i] = amps[i] - fit
    
  print('polifit baseline ' + str(i) + ' of ' + str(len(amps)))
  buffer = np.zeros(len(amps[0]))
  avg = np.median(amps[i])
  std = np.std(amps[i])
  buffer = (amps[i] > avg+3*std).astype(int)
  occup += buffer

occup = 100.0 / len(amps) * occup

print ("Plotting...")

fig, ax = plt.subplots()
ax.plot(freqs, occup, 'r', label='Occupancy', linewidth=1)
#ax.plot(freqs[-1], aspect, 'g', label='Average', linewidth=1)
ax.set_xlabel("Frequency [MHz]")
ax.set_ylabel("Occupancy [%]");
ax.set_autoscaley_on(False)
ax.set_ylim([0, 100.0])
ax.set_autoscalex_on(False)
ax.set_xlim([freqs[0], freqs[-1]])
ax.grid(1)
fig.set_size_inches(30, 20)

fgnm = 'rfispectrum_occupancy_full'
print ("Outputting spectrum occupancy plot %s" % fgnm)
fout = '%s/%s' % (outdir, fgnm)
fig.savefig(fout, dpi=600, bbox_inches='tight')
plt.close()

#plt.figure()
#plt.plot(freqs,np.transpose(amps))
#plt.title('Polyfit baseline all measurements')

#%% Occupancy in FDV's method:


print ("Calculating FDV occupancy...")
occup = np.zeros(len(D))

freqs = np.array(freq[(freq>=fmin) & (freq<=fmax)])
amps = np.array(D)

# calculate the envelope of the data:
mini = np.min(amps,0)
# remove the influence of RFI present 100% of the time:
#Frange=[360,380]
#A_aux = np.array( (mini[freqs==360],mini[freqs==365],mini[freqs==370],mini[freqs==375],mini[freqs==380] ))
#f_aux = np.array([360, 365, 370, 375, 380])
#A_interp = np.interp(freqs[(freqs>=360) & (freqs<=380)],f_aux,np.reshape(A_aux,5))
#mini[(freqs>=360) & (freqs<=380)] = A_interp


# remove the influence of RFI present 100% of the time:
#divide the frequency range in 5 MHz chunks
step =5
N = int((freqs[-1]-freqs[0])/step)
for i in range(N):
    Frange=[freqs[0]+i*step,freqs[0]+i*step+step]
#    Frange=[freqs[i],freqs[i]+step]
    A_aux = mini[(freqs>=Frange[0]) & (freqs<=Frange[1])]
    F_aux = freqs[(freqs>=Frange[0]) & (freqs<=Frange[1])]
#    A_minmax = np.sort(A_aux)[-1::-1]
    threshold = 0.1
    points = A_aux[np.where((A_aux-np.min(A_aux))<=threshold)]
    points[0] = A_aux[0]
    points[-1] = A_aux[-1]
    f_points = F_aux[np.where((A_aux-np.min(A_aux))<=threshold)]
    f_points[0] = F_aux[0]
    f_points[-1] = F_aux[-1]
    
    A_interp = np.interp(freqs[(freqs>=Frange[0]) & (freqs<=Frange[1])],f_points,points)
    mini[(freqs>=Frange[0]) & (freqs<=Frange[1])] = A_interp



amps_min = amps - mini

occup_avg = np.zeros(len(amps[0]))
occup_median = np.zeros(len(amps[0]))
occup_thresh = np.zeros(len(amps[0]))

for i in range(len(amps)):
#for k in range(20): #for debug
#    i = 200+k # for debuging
    thresh = mini - np.average(mini) + np.average(amps[i]) + np.std(amps_min[i])*1.5
    avg = np.average(amps_min[i])
    median = np.median(amps_min[i])
    std = np.std(amps_min[i])
    
    aux_avg = (amps_min[i] > avg+3*std).astype(int)
    aux_median = (amps_min[i] > median+3*std).astype(int)
    aux_thresh = (amps[i] > thresh).astype(int)

    occup_avg += aux_avg
    occup_median += aux_median
    occup_thresh += aux_thresh

#    plt.figure() #for debug
#    plt.plot(freqs,amps[i])
#    plt.plot(freqs,thresh)
#    plt.plot(freqs,mini)

    print('min vlue baseline ' + str(i) + ' of ' + str(len(amps)))


occup_avg = 100.0 / len(amps) * occup_avg
occup_median = 100.0 / len(amps) * occup_median
occup_thresh = 100.0 / len(amps) * occup_thresh

print ("Plotting...")

fig, ax = plt.subplots()
ax.plot(freqs, occup_avg, 'r', label='Occupancy over average', linewidth=1)
ax.plot(freqs, occup_median, 'b', label='Occupancy over median', linewidth=1)
ax.plot(freqs, occup_thresh, 'g', label='Occupancy over threshold', linewidth=1)
ax.set_xlabel("Frequency [MHz]")
ax.set_ylabel("Occupancy [%]");
ax.set_autoscaley_on(False)
ax.set_ylim([0, 100])
ax.set_autoscalex_on(False)
ax.set_xlim([freqs[0], freqs[-1]])
ax.grid(1)

fig.set_size_inches(30, 20)


fgnm = 'rfispectrum_occupancy_full_FDV'
print ("Outputting spectrum occupancy plot %s" % fgnm)
fout = '%s/%s' % (outdir, fgnm)
fig.savefig(fout, dpi=600, bbox_inches='tight')
plt.close()

#plt.figure()
#plt.plot(freqs,np.transpose(amps_min))
#plt.title('min baseline all measurements')


#plt.figure()
#plt.plot(freqs,amps[i])
#plt.plot(freqs,thresh)
#plt.plot(freqs,mini)


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



