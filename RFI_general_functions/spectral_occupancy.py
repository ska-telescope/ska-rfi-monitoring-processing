# -*- coding: utf-8 -*-
"""
Created on Wed May 29 12:32:37 2019
Occupancy in FDV's method:
# doing this in one hour chunks considering the variability of the background noise.
@author: F.Divruno
"""
import numpy as np
import matplotlib.pyplot as plt

def spectral_occupancy(freqs,D,outdir,title,std_multiplier):
        
    occup_thresh = np.zeros(len(D[0]))
    
    amps = np.array(D)
    tr = np.ndarray([np.size(amps,0),np.size(amps,1)])
    n_amps = len(amps)
    for j in range(n_amps):    
        blfit = 2
        A = amps[j]
        N = len(A)
        fit = np.zeros(N)
        K = int((freqs[-1]-freqs[0])/12)  #divide the frequency range in peaces of 12 MHz
        for z in range(K):
            freqs2 = freqs[z*int(N/K):(z+1)*int(N/K)]
            A2 = A[z*int(N/K):(z+1)*int(N/K)]
            coeffs = np.polyfit(freqs2,A2,blfit)
            fit[z*int(N/K):(z+1)*int(N/K)] = np.poly1d(coeffs)(freqs2)
        fit[-1]= fit[-2]
        err = (A - fit) # normalized spectrum
        tr = np.average(err)+np.std(err)*std_multiplier #threshold
        aux_thresh = (err > tr).astype(int)
        
        occup_thresh += aux_thresh
        
            # ----- debug ----
#            plt.figure() #for debug
#            plt.plot(freqs,amps[i])
#            plt.plot(freqs,thresh)
#            plt.plot(freqs,mini)
#            plt.legend(['amps','thresh','mini'])
#            plt.savefig(outdir+'debug_occup_'+str(k), dpi=100, bbox_inches='tight')
#
#            plt.figure() #for debug
#            plt.plot(freqs,10*np.log10(amps[i]))
#            plt.plot(freqs,10*np.log10(thresh))
#            plt.plot(freqs,10*np.log10(mini))
#            plt.legend(['amps','thresh','mini'])
#            plt.savefig(outdir+'debug_occup_log_'+str(k), dpi=100, bbox_inches='tight')
            # ----- debug ----
            
        print('min value baseline ' + str(j) + ' of ' + str(len(amps)))
    
    occup_thresh = 100.0 / len(D) * occup_thresh
    print ("Plotting...")
    
    fig, ax = plt.subplots()
    ax.plot(freqs, occup_thresh, 'g', label='Occupancy over threshold', linewidth=1)
    ax.set_xlabel("Frequency [MHz]")
    ax.set_ylabel("Occupancy [%]");
    ax.set_autoscaley_on(False)
    ax.set_ylim([0, 100])
    ax.set_autoscalex_on(False)
    ax.set_xlim([freqs[0], freqs[-1]])
    ax.grid(1)
    
    fig.set_size_inches(30, 20)
    
    
    fgnm = 'RFI_occupancy_%d to %d' % (int(freqs[0]),int(freqs[-1]))
    print ("Outputting spectrum occupancy plot %s" % fgnm)
    fout = '%s%s' % (outdir, fgnm)
    fig.savefig(fout, dpi=600, bbox_inches='tight')
#    plt.close()
    
    #plt.figure()
    #plt.plot(freqs,np.transpose(amps_min))
    #plt.title('min baseline all measurements')
    
    
    #plt.figure()
    #plt.plot(freqs,amps[i])
    #plt.plot(freqs,thresh)
    #plt.plot(freqs,mini)
    return occup_thresh
