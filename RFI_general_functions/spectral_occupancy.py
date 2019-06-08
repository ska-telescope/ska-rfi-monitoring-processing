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
    f_step = freqs[1]-freqs[0]
    BW = freqs[-1]-freqs[0]
    N = len(freqs)
    
    #find parts of the spectrum with significant power
    Ave = np.average(10**(D/10),0)
    channel = 1 #1 MHz channel
    points_in_channel = int(channel/f_step)
    channels_in_BW = int(N/points_in_channel)
    Pow_in_ch = np.zeros(channels_in_BW)
    freq_in_ch = np.zeros(channels_in_BW)
    for i in range(channels_in_BW):
        Pow_in_ch[i] = 10*np.log10(np.sum(Ave[i*points_in_channel:(i+1)*points_in_channel]))
        freq_in_ch[i] = freqs[i*points_in_channel]
    dif = Pow_in_ch[1::]-Pow_in_ch[0:-1:1]
    dif = np.append(dif,0)
#    dif = np.append(dif,0)
    f_RFI_start = [] #np.ndarray([])
    f_RFI_stop = [] #np.ndarray([])
    rfi_flag = 0
    for i in range(len(dif)-1):
        if (dif[i]>=1) & (rfi_flag==0):
            if i>0:
                f_RFI_start = np.append(f_RFI_start,freq_in_ch[i-1])
                rfi_flag = 1
        if (dif[i]<=-1) & (rfi_flag==1):
            if i<len(dif):
                rfi_flag = 2
        if rfi_flag == 2:
            if dif[i]>-1:
                f_RFI_stop = np.append(f_RFI_stop,freq_in_ch[i])                
                rfi_flag = 0
    
#    plt.plot(freq_in_ch,dif)            
#    plt.stem(f_RFI_start,np.ones(len(f_RFI_start)))
#    plt.stem(f_RFI_stop,np.ones(len(f_RFI_stop))*-1)
    ind_start = np.zeros(len(f_RFI_start)).astype('int')
    ind_stop = np.zeros(len(f_RFI_start)).astype('int')
    for i in range(len(f_RFI_start)):
        ind_start[i] = int(np.where(freqs>=f_RFI_start[i])[0][0])
        ind_stop[i] = int(np.where(freqs<=f_RFI_stop[i])[0][-1])
        
    
    occup_thresh = np.zeros(len(D[0]))
    
    amps = np.array(D)
    err = np.ndarray([np.size(amps,0),np.size(amps,1)])
    n_amps = len(amps)
    for j in range(n_amps):  
        blfit = 2
        A = amps[j]
        A_aux = np.array(A)
        for i in range(len(f_RFI_start)): #gets rid of the RFI present all the time for doing the fit.
#            A_aux[(freqs>=f_RFI_start[i])&(freqs<=f_RFI_stop[i])] = np.ones(sum((freqs>=f_RFI_start[i])&(freqs<=f_RFI_stop[i])))*(A[(freqs>=f_RFI_start[i])][0]+A[(freqs<=f_RFI_stop[i])][-1])/2 
            A_aux[ind_start[i]:ind_stop[i]] = np.ones(ind_stop[i]-ind_start[i])*(A[ind_start[i]]+A[ind_stop[i]])/2 
            
        fit = np.zeros(N)
        #K = int(/15)  #divide the frequency range in peaces of 12 MHz
        N_12 = int(12/f_step) # aprox number of steps to get 12 MHz
        index = int(0)
        while (index+N_12)<N:
            freqs2 = freqs[index:index+N_12]
            A2 = A_aux[index:index+N_12]
            coeffs = np.polyfit(freqs2,A2,blfit)
            fit[index:index+N_12] = np.poly1d(coeffs)(freqs2)
            index += N_12
#        in the case there are more points
        freqs2 = freqs[index::]
        A2 = A_aux[index::]
        coeffs = np.polyfit(freqs2,A2,blfit)
        fit[index::] = np.poly1d(coeffs)(freqs2)
        
        err[j] = (A - fit) # normalized spectrum
        tr = np.average(err[j])+np.std(err[j])*std_multiplier #threshold
        aux_thresh = (err[j] > tr).astype(int)
        
        occup_thresh += aux_thresh
        
            # ----- debug ----
        if j ==1:
            plt.figure() #for debug
            plt.plot(freqs,A)
            plt.plot(freqs,fit)
            plt.legend(['amps','fit'])
            plt.savefig(outdir+'debug_occup_'+str(j), dpi=100, bbox_inches='tight')
            
            
            plt.figure() #for debug
            plt.plot(freqs,err[j])
            plt.plot(freqs,tr*np.ones(len(freqs)))
            plt.legend(['error','thresh'])
            plt.savefig(outdir+'debug_occup_'+str(j), dpi=100, bbox_inches='tight')

            plt.figure() #for debug
            plt.plot(freqs,10*np.log10(err[j]))
            plt.plot(freqs,10*np.log10(tr*np.ones(len(freqs))))
            plt.legend(['amps','thresh'])
            plt.savefig(outdir+'debug_occup_log_'+str(j), dpi=100, bbox_inches='tight')
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
    return err,occup_thresh
