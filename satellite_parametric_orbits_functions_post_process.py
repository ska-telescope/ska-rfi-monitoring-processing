# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 01:26:31 2019

@author: f.divruno
"""


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os as os
import astropy.units as u

#%% Functions

def import_data(directory,planes,sats,idx=0):
    files = os.listdir(directory)
#    find the files
    valid_files = list([])
    for file in files:
        if file[-4::]== '.npz':
            N_planes = float(file[file.find(' planes')-2:file.find(' planes')])
            Sats_per_plane = float(file[file.find(' sats pp')-2:file.find(' sats pp')])
            if (N_planes == planes)  & (Sats_per_plane==sats):
                valid_files.append(file)
#    Load the values                
    if len(valid_files)>0:
        aux = np.load(directory+'/'+valid_files[idx])
        el = aux['el']
        az = aux['az']
        seed = int(valid_files[idx][valid_files[idx].find(' seed ')+6:valid_files[idx].find(' seed ')+10])
        correct_dBm = 1 #some calculations with error in dBm to dBW, needs to add 30 dB
        if correct_dBm:
            Prx = aux['Prx']*1e3
            maxPrx = aux['maxPrx']*1e3
            avePrx = aux['avePrx']*1e3
        else:
            Prx = aux['Prx']
            maxPrx = aux['maxPrx']
            avePrx = aux['avePrx']
        print('\n Loaded file: "%s"'%(valid_files[idx]))
        return(seed,el,az,Prx,maxPrx,avePrx,len(valid_files))

    else:
        print('\nERROR:\n Specified number of planes or sats not found')
        return(0,0,0,0,0,0,0)

def plot_trimesh(el,az,Z,title='', view_el = 60, view_az = -90):
    '''
    
    '''
    
    import matplotlib.tri as mtri                
    y = np.array(el)
    x = np.array(az)
    z = 10*np.log10(np.array(Z))
     
    triang = mtri.Triangulation(x, y)
    
    fig = plt.figure(figsize=[15,9])
    ax = fig.add_subplot(1,1,1,projection='3d')
    
    ax.plot_trisurf(triang, z, cmap='jet', edgecolor='none', linewidth=0, antialiased=False)
    
    ax.view_init(elev = view_el, azim = view_az)
    
    ax.set_xlabel('Azimuth')
    ax.set_ylabel('Elevation')
    ax.set_zlabel('Averaged power dBm ')
    plt.title(title)
    return fig,ax
    

def plot_rx_power_in_time(Point_az, Point_el, Prx_time,fig=[]):
    
    ind = np.where((el>=Point_el) & (az>= Point_az))[0][0]
    
    if fig==[]:
        plt.figure()
    else:
        plt.figure(fig)
    plt.plot(10*np.log10(Prx_time[ind]+1e-20))
    plt.title('received power in time poining: el= %.1f deg, az = %f deg'%(Point_el, Point_az))
#    plt.title('time domain received power in Az:%f , El: %f' %(Point_az,Point_el))


#%% Start rhe calculation
if __name__ == '__main__':
    
    N_planes = 24
    Sats_per_plane = 66 
    
    [seed,el,az,Prx,maxPrx,avePrx,N_files] = import_data('../satellite_results',N_planes,Sats_per_plane)
    
    seed_v = np.zeros(N_files)
    Prx_v = np.zeros([len(Prx[:,0]),len(Prx[0,:]),N_files])
    maxPrx_v = np.zeros([len(maxPrx),N_files])
    avePrx_v = np.zeros([len(avePrx),N_files])

#    N_files = 1    
    for i in range(N_files):
                
        max_time = len(Prx[0,:])
        
        identifier = '- %d planes - %d sats pp - seed %d '%(N_planes,Sats_per_plane,seed) # for plotting and saving
        
        
        # plot the received power in the sky
        blk,ax = plot_trimesh(el,az,maxPrx,'Maximum received power '+ identifier ,0,0)
        plt.savefig('../satellite_results/Max received power - full sky '+identifier+'-side.png')
        ax.view_init(90,-90)
        plt.draw()
        plt.savefig('../satellite_results/Max received power - full sky '+identifier+'-front.png')
        plt.close()
                
        blk,ax = plot_trimesh(el,az,avePrx, 'Average received power in %s'%(max_time)+ identifier ,0,0)
        plt.savefig('../satellite_results/Avg received power - full sky '+identifier+'-side.png')
        ax.view_init(90,-90)
        plt.draw()
        plt.savefig('../satellite_results/Avg received power - full sky '+identifier+'-front.png')
        plt.close()
        
        #max in time domain:
        k = np.where(maxPrx==np.max(maxPrx))
        plot_rx_power_in_time(az[k],el[k],Prx)
        plt.savefig('../satellite_results/Instantaneous received power - el %.2f Az %.2f'%(el[k],az[k])+identifier+'.png')
        plt.close()
        
        seed_v[i] = seed
        Prx_v[:,:,i] = Prx
        maxPrx_v[:,i] = maxPrx
        avePrx_v[:,i] = avePrx
        
        if i < N_files-1:
            [seed,el,az,Prx,maxPrx,avePrx,blk] = import_data('../satellite_results',N_planes,Sats_per_plane,idx=i+1)

#%%    Clculate metrics in function of elevation
    el_values = np.unique(el)
    max_maxPrx = np.zeros([len(el_values),(N_files)])
    ave_maxPrx = np.zeros([len(el_values),(N_files)])
    max_avePrx = np.zeros([len(el_values),(N_files)])
    ave_avePrx = np.zeros([len(el_values),(N_files)])

    for i in range(N_files):
        for k in range(len(el_values)):
            ind = np.where(el==el_values[k])
            max_maxPrx[k,i] = np.max(maxPrx[ind])
            ave_maxPrx[k,i] = np.mean(maxPrx[ind])
            max_avePrx[k,i] = np.max(avePrx[ind])
            ave_avePrx[k,i] = np.mean(avePrx[ind])
#%% 
    identifier2 = '- %d planes - %d sats pp '%(N_planes,Sats_per_plane) # for plotting and saving, the random seed doesnt matter here        
    plt.figure()
    plt.plot(el_values,10*np.log10(max_maxPrx))
    plt.title('Max values in elevation of MaxPrx')
    plt.xlabel('Elevation')
    plt.ylabel('dBm')    
    plt.savefig('../satellite_results/Max values in elevation of MaxPrx'+identifier2+'.png')
       
    plt.figure()
    plt.plot(el_values,10*np.log10(ave_maxPrx))
    plt.title('Average values in elevation of MaxPrx')
    plt.xlabel('Elevation')
    plt.ylabel('dBm')    
    plt.savefig('../satellite_results/Average values in elevation of MaxPrx'+identifier2+'.png')

    plt.figure()
    plt.plot(el_values,10*np.log10(max_avePrx))
    plt.title('Max values in elevation of AvePrx')
    plt.xlabel('Elevation')
    plt.ylabel('dBm')    
    plt.savefig('../satellite_results/Max values in elevation of AvePrx'+identifier2+'.png')
    
    plt.figure()
    plt.plot(el_values,10*np.log10(ave_avePrx))
    plt.title('Average values in elevation of AvePrx')
    plt.xlabel('Elevation')
    plt.ylabel('dBm')    
    plt.savefig('../satellite_results/Average values in elevation of AvePrx'+identifier2+'.png')    
        
#%% Calculate the number of hits 
    threshold = 10**(-100/10) # threshold in linear units [mW]
    N_hits = np.zeros([len(el_values),(N_files)])
    for i in range(N_files):
        for k in range(len(el_values)):
            ind = np.where(el==el_values[k]) # for one elevation value
            ind2 = np.where(Prx[ind,:].flatten()>=threshold) #all time samples at those elevation values
            N_hits[k,i] = len(ind2[0])/len(Prx[ind,:].flatten())*100
            
    plt.figure()
    plt.plot(el_values,N_hits)
    plt.title('number of hits (prx>%.1f) in elevation in simulation time'%(10*np.log10(threshold)))
    plt.xlabel('Elevation')
    plt.ylabel('Hits number in %')    
    plt.savefig('../satellite_results/Percentage of hits in elevation of Prx'+identifier2+'.png')
        
   
#%% calculate the histograms
    plt.figure(figsize=[15,10])
    cdf = plt.hist(10*np.log10(Prx.flatten()),1000,density=True,cumulative=True,histtype='step')
    plt.title('Cumulative dist funct of instantaneous received power')
    plt.savefig('../satellite_results/Cumulative dist funct of Prx'+identifier2+'.png')
    i = np.where(cdf[0]>0.98)[0][0]
    plt.plot([cdf[1][i],cdf[1][i]],[0,1])
    
    plt.figure(figsize=[15,10])
    cdf = plt.hist(10*np.log10(avePrx.flatten()),1000,density=True,cumulative=True,histtype='step')
    plt.title('Cumulative dist funct of average received power')
    plt.savefig('../satellite_results/Cumulative dist funct of avePrx'+identifier2+'.png')
    i = np.where(cdf[0]>0.98)[0][0]
    plt.plot([cdf[1][i],cdf[1][i]],[0,1])

    plt.figure(figsize=[15,10])
    cdf = plt.hist(10*np.log10(maxPrx.flatten()),1000,density=True,cumulative=True,histtype='step')
    plt.title('Cumulative dist funct of maximum received power')
    plt.savefig('../satellite_results/Cumulative dist funct of maxPrx'+identifier2+'.png')
    i = np.where(cdf[0]>0.98)[0][0]
    plt.plot([cdf[1][i],cdf[1][i]],[0,1])
