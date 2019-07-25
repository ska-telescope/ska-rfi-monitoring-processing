# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 22:15:32 2019

    Assesment of satellite constellation on SKA

@author: f.divruno
"""



from astropy import units as u
from astropy.time import Time
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import rfiLib as RFI
from rfiLib.Pointing_DISH import *

from mpl_toolkits.mplot3d import Axes3D
import cartopy.crs as ccrs
import pyproj


font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 24}

matplotlib.rc('font', **font)


Re = Earth.R
GHz = 1e9
MHz = 1e6
kHz = 1e3
km = 1e3

def rotation_matrix(axis, theta):
    '''
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    '''
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def Generate_NGSO(height,inc,raan,nu,rand_seed=0):
    '''
        Generate the orbit elements given the clasic parameters
    '''

    R_e = Earth.R
    a = R_e + height*u.km
    ecc = 0. * u.one
    inc = inc * u.deg
    raan = raan * u.deg
    argp = 0 * u.deg
    nu = nu * u.deg
    epoch = Time("2015-01-01 00:00:00", scale="tdb")
    
    np.random.seed(rand_seed)
    epoch += np.random.random(1)*3600
    
    return Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu,epoch)

def antenna_gain(angle):
    '''
        SKA antenna gain considering that is equal in azimuth, only varying in 
        elevation.
    '''
    
    N = np.size(angle)
    G = np.zeros(N)
    delay_samples = np.zeros(N)
    angle = np.abs(angle)
    if N >1:
        for i in range(N):
            if angle[i]<0.1:
                G[i] = 60
#                delay_samples[i] = 0
        
            elif (angle[i]>0.1) & (angle[i]<3.5):
                G[i] = 32-25*np.log10(angle[i]) # dB
#                delay_samples[i] = np.random.rand(1)*10 # maximum random delay because of the sidelobes = 10/fs = 2.5 ns en 4Gsps
                    
            else:
                G[i]=0
#                delay_samples[i] = np.random.rand(1)*10
    else:
        if angle<0.1:
            G = 60
#            delay_samples = 0
    
        elif (angle>0.1) & (angle<3.5):
            G = 32-25*np.log10(angle) # dB
#            delay_samples = np.random.rand(1)*10 # maximum random delay because of the sidelobes = 10/fs = 2.5 ns en 4Gsps
                
        else:
            G=0
#            delay_samples = np.random.rand(1)*10
 
    return delay_samples,G # gain in the beam    
    
def Free_space(d,f):
    '''
        Free Space Path Loss
    '''
    try:
        test1 = (d.unit=='m') & (f.unit=='MHz')
        if test1:
            FSPL = 20*np.log10(d.value) + 20*np.log10(f.value) - 25.77 #height in meters, freq in MHz
        else:
            raise Exception

    except:
        print('\nERROR Free_space:\nNeed to use m and MHz units\n')
        raise SystemExit
        
    return FSPL


ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')



def received_power(Rx_vect,T,Pointing_el,Pointing_az,A):
    '''
    Calculates the received power by a station in the position Rx_vect,
    Pointing to the direction el,az.
    Rx_vect:
    T: Transformation matrix from receptor ENU to ECEF 
    Pointing_el: elevation pointing
    Pointing_az: azimuth pointing
    A: position of each satellite in [sat number,time step, [x,y,z]]
    '''
    
    
    
    N_sats = np.size(A,0)
    freq = 11*u.GHz
    
    pointing = [Pointing_el,Pointing_az]*u.deg #  elevation, Azumith
    
    
    
    # transform the pointing direction to the ECEF frame
    pointingENU = np.array(Pointing_to_ENU(pointing[0],pointing[1]))
    pointingECEF = np.matmul(T,pointingENU) #unit vector

    
    PSDrx = np.ones([np.size(A,1),N_sats])*-500 
      
    for i in range(N_steps):
        D = np.transpose(A[:,i,:]) - Rx_vect
        D_norm = np.linalg.norm(D,axis=0)
        index = np.where(D_norm<2700)[0].astype('int32')
        
        FSPL = Free_space((D_norm[index]*u.km).to('m'),freq.to('MHz'))
        
        Point_vect = np.dot(np.reshape(np.array(pointingECEF),[3,1]),np.ones([1,len(index)]))
        cos_alpha1 = ((Point_vect*D[:,index]).sum(0)/D_norm[index])    
        
        alpha = np.arccos(cos_alpha1)*180/np.pi
       
        [blk, G_SKA] = antenna_gain(alpha)
   
        PSDrx[i,index] = EIRP_Hz + 30 - FSPL + G_SKA + G_sat_sidelobes #30 is to pass EIRP to dBm/Hz 

    
    PSDrx_t = 10*np.log10(np.sum(10**(PSDrx/10),1)) # sum of all the satellites contributions for each timestep
        
    return PSDrx_t  


def propagate_orbits(step,N_steps,sats,rand_seed=0):
    '''
        step: time step
        N_steps: number of time steps
        sats: list of satellite orbits
        
    '''

    
    N_sats = N_planes*Sats_per_plane
    angle_vel = 2*np.pi/24/3600 #rad/sec
    theta = 0
    axis =[0,0,1] # Earth obliquity is 23deg, should this be considered?    

    np.random.seed(rand_seed)
    rand_start = np.random.rand(1)*3600*u.s
    
    A = np.zeros([N_sats,N_steps,3])
    i=0
    for sat in sats:

        for j in range(N_steps):
            theta += step*angle_vel    
            ss_delta = sat.propagate(step*u.s + rand_start) #propagate the orbit for the step time
            [x,y,z] = ss_delta.r # get the position in ECEF frame
            [x1,y1,z1] = np.dot(rotation_matrix(axis, theta), [x.value,y.value,z.value]) # Rotate the frame along Z axis (simulate earth rotation)
            A[i,j,:] = np.array([x1,y1,z1])
            print('Sat number %d - Timestep %d' %(i,j))    
        i+=1
    return A

#%% satellite parameters

N_planes = 24
Sats_per_plane = 66
N_sats = int(N_planes*Sats_per_plane)
raan_step = 360/N_planes        #Right ascencion steps
nu_step = 360/Sats_per_plane    #Declination steps
height = 550    
inclination = 53    


G_sat_sidelobes = -10 # sidelobes

P_saturation = -90 #dBm
R_earth = 6300*km


# satellite parameters
EIRP_4kHz = -15 #dBW/4kHz
BW_ch = 250*MHz
fmin = 10.7*GHz #minimum assigned frequency
fmax = 12.7*GHz #max assigned frequency


#calcs
freq = 11*u.GHz # np.linspace(fmin,fmax,20)
EIRP_Hz = EIRP_4kHz + 10*np.log10(1/(4*kHz))
EIRP_250MHz = EIRP_Hz + 10*np.log10(250*MHz)
EIRP_2000MHz = EIRP_Hz + 10*np.log10(2000*MHz)



#%% Simulation control: 
# Core location
SKA_core_lat = -30.71278*u.deg
SKA_core_lon = 21.44305*u.deg
SKA_core_h = 1000*u.m
SKA_core_ECEF = np.array(pyproj.transform(lla,ecef,SKA_core_lon.value,SKA_core_lat.value,SKA_core_h.value))/km

#Pointing parameters
elev_min = 15*u.deg
elev_max = 90*u.deg
Az_min = 15*u.deg
Az_max = 90*u.deg

#Time control
step = 30  # time step for the propagation of the satellites.
totalHours = 0.1 #total duration of the simulation


#%% Generate the orbits for each satellite.

totalTime = totalHours*3600
N_steps = int(totalTime/step)

try: #Try to load a saved file.
    Aux = np.load('Starlink_constellation_propagation_%dhrs_step%d.npz'%(totalHours,step))
    A = Aux['A']
    step = Aux['step']
    raise Exception
except:
    
    sats = list()
    for i in range(N_planes):
        for k in range(Sats_per_plane):
            sats.append(Generate_NGSO(height,inclination,i*raan_step,k*nu_step))
            print(i,k)
    
    A = propagate_orbits(step,N_steps,sats,rand_seed=0)    
    
#    np.savez('Starlink_constellation_propagation_%dhrs_step%d'%(totalHours,step),A=A,step=step)     

time = np.linspace(0,step*(N_steps-1),N_steps)

#%% Plot orbits in the full time calculated
   # plot over the map with PLateeCarree projection.

plot_Map = 0
if plot_Map:    
    fig = plt.figure(figsize=[20,10])
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.stock_img()
    
    for i in range(N_steps):
        lon,lat,alt = pyproj.transform(ecef,lla,A[:,i,0]*1e3,A[:,i,1]*1e3,A[:,i,2]*1e3)
        ax.plot(lon,lat, color='blue', linewidth=0.1,transform=ccrs.Geodetic())
        print('Drawing trajectory, completed: %.2f'%(i/N_steps*100))



#%% Calculate the agregated power as a function of time
# Consier case of SKA-MID pointing to zenith and  the NGSO not pointing to SKA antenna


pointing = [45,90]*u.deg #  elevation, Azumith


# transform the pointing direction to the ECEF frame
pointingENU = np.array(Pointing_to_ENU(pointing[0],pointing[1]))
T = ENU_to_ECEF(SKA_core_lat,SKA_core_lon) 
pointingECEF = np.matmul(T,pointingENU) #unit vector



pointing_lat,pointing_lon = ENU_to_latlon(pointingECEF*(Re.value+550e3))
# PLot the distribution of the satelites with the GAin as colors
fig = plt.figure(figsize=[15,8])
ax = plt.axes(projection=ccrs.PlateCarree())
ax.stock_img()

#SKA_core_lon,SKA_core_lat,alt = pyproj.transform(ecef,lla,SKA_core_ECEF[0]*1e3,SKA_core_ECEF[1]*1e3,SKA_core_ECEF[2]*1e3)
ax.plot(SKA_core_lon,SKA_core_lat,linewidth=4,marker='X',transform=ccrs.Geodetic())
ax.plot(pointing_lon,pointing_lat,linewidth=4,marker='d',transform=ccrs.Geodetic())

## Initializing For

PSDrx = np.ones([len(time),int(N_sats)])*-500 
alphaList = list()
SKA_core_norm = np.linalg.norm(SKA_core_ECEF)
Sat_norm = np.linalg.norm(A[0,0,:]) #all satellites are the same distance from the centre of the earth.
SKA_core_vect = np.dot(np.reshape(np.array(SKA_core_ECEF),[3,1]),np.ones([1,N_sats]))

# ==================== DEBUG
#N_steps = 1  # force the time steps to one.

# ==================== DEBUG
    
g_list1 = list()
for i in range(N_steps):
    D = np.transpose(A[:,i,:])-SKA_core_vect
    D_norm = np.linalg.norm(D,axis=0)*u.km
    index = np.where(D_norm<2700*u.km)[0].astype('int32')
    
    FSPL = Free_space(D_norm[index].to('m'),freq.to('MHz'))
    
    Point_vect = np.dot(np.reshape(np.array(pointingECEF),[3,1]),np.ones([1,len(index)]))
    cos_alpha1 = ((Point_vect*D[:,index]).sum(0)/D_norm[index].value)    

    alpha = np.arccos(cos_alpha1)*180/np.pi
    
    alphaList.append(alpha)
    [blk, G_SKA] = antenna_gain(alpha)

    satLon1,satLat1,alt = pyproj.transform(ecef,lla,A[index,i,0]*1e3,A[index,i,1]*1e3,A[index,i,2]*1e3)
    ax.scatter(satLon1,satLat1,s=10,c = 10**(G_SKA/10),transform=ccrs.Geodetic(),cmap='jet',alpha=1)    
#    ax.scatter(satLon1,satLat1,s=10,c =(1-alpha/90),transform=ccrs.Geodetic(),cmap='jet',alpha=1)    

    g_list1.append(G_SKA)
    print('Calculating received power, completed: %.2f'%(i/N_steps*100))
    PSDrx[i,index] = EIRP_Hz + 30 - FSPL + G_SKA + G_sat_sidelobes #30 is to pass EIRP to dBm/Hz 

plt.title('Angle distribution- pointing: elev= %.2f  Az= %.2f '%(pointing[0].value,pointing[1].value))

PSDrx_t = 10*np.log10(np.sum(10**(PSDrx/10),1)) # sum of all the satellites contributions for each timestep
    

#SKA_std = RFI.SKA_EMIEMC_std(freq,freq,20,0)

plt.figure()
plt.plot(time,PSDrx_t)
plt.title('Equivalent PSD received every %d seconds'%(step))
plt.ylabel('dBm/Hz')

plt.figure()
plt.plot(time,PSDrx_t + 10*np.log10(BW_ch))
plt.title('Power received every %d seconds in 250 MHz BW'%(step))
plt.ylabel('dBm')
plt.plot([time[0],time[-1]],[-90,-90])


Prx = PSDrx_t + 10*np.log10(250*MHz)
count = np.sum(np.array(Prx > (P_saturation)).astype('int'))
time_perc = count/len(Prx)*100
print('Percentage of the time the power is above the integrated noise power: %d %%'%(time_perc))


#%% Creating a grid in the sky:

SKA_core_vect = np.dot(np.reshape(np.array(SKA_core_ECEF),[3,1]),np.ones([1,N_sats]))
T = ENU_to_ECEF(SKA_core_lat,SKA_core_lon) 

elev_step = 10
elev = np.linspace(elev_step,elev_max.value,int((elev_max.value-elev_min.value)/elev_step))*u.deg
N_elev = len(elev)
PSD_rx = list()
el = list()
az = list()


for i in range(len(elev)):
    width = np.array(90/N_elev/np.cos(elev[i]*np.pi/180)).astype('int32')
    Az = (np.linspace(-180,180,int(360/width)))*u.deg
   
    for j in range(len(Az)):
        psd_rx = received_power(SKA_core_vect,T,elev[i],Az[j],A)
        PSD_rx.append(psd_rx)
        el.append(elev[i])
        az.append(Az[j])
        print('Calculating El: %s  - Az: %s'%(elev[i],Az[j]))
    




#%%        
#generate the grid to plot
el2 = np.zeros(len(el))
az2 = np.zeros(len(az))
for i in range(len(el)):
    el2[i] = el[i].value        
    
for i in range(len(az)):
    az2[i] = az[i].value        
    
y = np.array(el2)
x = np.array(az2)
z = np.array(PSD_rx+10*np.log10(BW_ch))
 

import matplotlib.tri as mtri

triang = mtri.Triangulation(x, y)

fig = plt.figure()
ax = fig.add_subplot(1,1,1,projection='3d')

ax.plot_trisurf(triang, z, cmap='jet')

ax.view_init(elev=90, azim=-90)

ax.set_xlabel('Azimuth')
ax.set_ylabel('Elevation')
ax.set_zlabel('Averaged power dBm ')

plt.show()

