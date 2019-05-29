# -*- coding: utf-8 -*-
"""
Created on Wed May 29 12:32:37 2019
Occupancy in FDV's method:
# doing this in one hour chunks considering the variability of the background noise.
@author: F.Divruno
"""
import numpy as np
import matplotlib.pyplot as plt

def spectral_occupancy(freqs,D,outdir,std_multiplier):
    N_chunk = int(len(D)/24)
    print ("Calculating occupancy...")
    
    occup_thresh = np.zeros(len(D[0]))
    
    for i in range(25):
        try:
            D1 = D[i*N_chunk:(i+1)*N_chunk]        
            print ("Chunk : " + str(i) + ' of 24')     
        except:
            D1 = D[i*N_chunk::]
            print ("Last chunk")     
           
        
        amps = np.array(D1)
        
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
        
        
        for i in range(len(amps)):
        #for k in range(20): #for debug
        #    i = 200+k # for debuging
            thresh = mini - np.average(mini) + np.average(amps[i]) + np.std(amps_min[i])*std_multiplier
#            avg = np.average(amps_min[i])
#            median = np.median(amps_min[i])
#            std = np.std(amps_min[i])
            
#            aux_avg = (amps_min[i] > avg+3*std).astype(int)
#            aux_median = (amps_min[i] > median+3*std).astype(int)
            aux_thresh = (amps[i] > thresh).astype(int)
        
#            occup_avg += aux_avg
#            occup_median += aux_median
            occup_thresh += aux_thresh
        
        #    plt.figure() #for debug
        #    plt.plot(freqs,amps[i])
        #    plt.plot(freqs,thresh)
        #    plt.plot(freqs,mini)
        
            print('min value baseline ' + str(i) + ' of ' + str(len(amps)))
    
    
#    occup_avg = 100.0 / len(D) * occup_avg
#    occup_median = 100.0 / len(D) * occup_median
    occup_thresh = 100.0 / len(D) * occup_thresh
    
    print ("Plotting...")
    
    fig, ax = plt.subplots()
#    ax.plot(freqs, occup_avg, 'r', label='Occupancy over average', linewidth=1)
#    ax.plot(freqs, occup_median, 'b', label='Occupancy over median', linewidth=1)
    ax.plot(freqs, occup_thresh, 'g', label='Occupancy over threshold', linewidth=2)
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
    fig.savefig(fout, dpi=300, bbox_inches='tight')
    plt.close()
    
    #plt.figure()
    #plt.plot(freqs,np.transpose(amps_min))
    #plt.title('min baseline all measurements')
    
    
    #plt.figure()
    #plt.plot(freqs,amps[i])
    #plt.plot(freqs,thresh)
    #plt.plot(freqs,mini)
    return occup_thresh