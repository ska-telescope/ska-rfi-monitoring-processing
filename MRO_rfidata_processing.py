# -*- coding: utf-8 -*-
"""
Created on Wed May 22 02:42:48 2019

@author: F.Divruno
"""
import matplotlib.pyplot as plt
import RFI_general_functions as RFI
import numpy as np

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype='float32')
    ret[n:] = ret[n:] - ret[:-n]
    ret[0:n-1] = a[0:n-1]*n #try
    #return ret[n - 1:] / n
    return ret/ n


#folderpath = r'C:\Users\F.Divruno\Dropbox (SKA)\14- RFI environment\01- Australia\rfidata_mro'
folderpath = r'/mnt/data/Downloads/datafiles'
filepath = folderpath + r'\rfidata_mro_2019-05-12_2340.fits.gz'

read_files = 0
if read_files:
    [f,a] = RFI.read_MRO_data(folderpath)
    np.savez_compressed(folderpath + r'/MRO_rfidata', f=f, A=a)
else:
    Aux = np.load(folderpath+r'\MRO_rfidata.npz')
    a = Aux['A']/10-107 #in V**2 originally, scaled to get to dBm
    f = Aux['f'] # in MHz   

#%%
n_ave = 20
b = np.array(a)
d = np.array(a)
c = np.array(a)
#c = np.zeros([np.size(a,0),np.size(a,1)-(n_ave-1)])

for i in range(len(a)):
    blfit = 3
    if blfit != 0:
        # subtract an Nth order polynomial to remove baseline fluctuations
        coeffs = np.polyfit(f,a[i],blfit)
        fit = np.poly1d(coeffs)(f)
#        print ("max spect residual: %0.1f" % np.sqrt(np.median((a[i]-fit)**2)))
        b[i] = a[i] - fit
        print (i)

#
#    ave = moving_average(a[i,:],n_ave)
##    c[i,:] = a[i,19:] - ave
#    c[i,:] = a[i,:] - ave
    
    miin = np.min(a,0)
    d[i] = a[i]- miin
    

#%%
    
        
plt.figure()
plt.plot(f,np.transpose(a[-1,:]))
plt.plot(f,fit)
plt.xlabel('Freq')
plt.ylabel('power dBm')

#plt.figure()
#plt.plot(np.transpose(a[-1,:]))
#plt.plot(ave)


        
plt.figure()
plt.plot(f,np.transpose(b[1,:]))
plt.plot(f,np.transpose(d[1,:]))
#plt.plot(f_ave,np.transpose(c[1,:]))
#plt.plot(f,np.transpose(c[1,:]))

#%%
ave = np.average(b,0)
plt.figure()
#plt.plot(f_ave,ave)
plt.plot(f,ave)

#%% 
RFI.plot_percentile(f,d,100,'dBm','MRO data, normalized with the minimum [dB]')
RFI.plot_percentile(f,b,100,'dBm','MRO data, normalized with polyfit [dB]')
RFI.plot_percentile(f,a,100,'dBm','MRO data [dBm]')
