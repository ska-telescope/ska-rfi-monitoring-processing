# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 01:26:31 2019

@author: f.divruno
"""


import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import cartopy.crs as ccrs
import astropy.constants as const


#%% Functions 
def antenna_gain_times(angle):
    '''
        SKA antenna gain considering that is equal in azimuth, only varying in 
        elevation.
    '''
 
    N = len(angle)
    G = np.zeros(N)

    angle = np.abs(angle)
    G[angle<=0.1] = 10**6
    G[(angle>0.1) & (angle<=3.5)] = 1.58e3*(angle[(angle>0.1) & (angle<3.5)]**(-2.5))
    G[angle> 3.5] = 1

    return G # gain in the beam   


def antenna_gain(angle):
    '''
        SKA antenna gain considering that is equal in azimuth, only varying in 
        elevation.
    '''
#    try:
#        if angle.unit == u.deg:
#            angle = angle.value
#        elif angle.unit == u.rad:
#            angle = angle.to(u.deg).value
#    except:
#        angle = angle
    
    dim0 = 0
    dim1 = 0
    try:
        dim0 = np.size(angle,0)
        dim1 = np.size(angle,1)
        angle = angle.flatten()
    except:
        angle = angle
        
    N = np.size(angle)
    G = np.zeros(N)
#    delay_samples = np.zeros(N)
    angle = np.abs(angle)
    if N >1:
        G[angle<0.1] = 60
        G[(angle>0.1) & (angle<3.5)] = 32-25*np.log10(angle[(angle>0.1) & (angle<3.5)])
        G[angle> 3.5] = 0
        
#        for i in range(N):
#            if angle[i]<0.1:
#                G[i] = 60
#        
#            elif (angle[i]>0.1) & (angle[i]<3.5):
#                G[i] = 32-25*np.log10(angle[i]) # dB
#                    
#            else:
#                G[i]=0
    else:
        if angle<0.1:
            G = 60
    
        elif (angle>0.1) & (angle<3.5):
            G = 32-25*np.log10(angle) # dB

                
        else:
            G=0
    
    if dim1>0:
        G = np.reshape(G,[dim0,dim1])
    return G # gain in the beam   

#%% create the orbits in an inertial frame as the GCRS
    
# define a circle of radius Earth.R + h in the GCRS ref system 

period = 96*u.min 
steps = 10000
N_orbits = 1
t = np.linspace(0,1,steps)*period.to(u.s)*N_orbits
epoch0 = Time("2015-01-01 00:00:00", scale="utc")
epoch = epoch0 + t
phi_max = 2*np.pi*N_orbits

N_planes = 24
Sats_per_plane = 66
N_sats = int(N_planes*Sats_per_plane)
raa_step = 360/N_planes*u.deg        #Right ascencion steps
mA_step = 360/Sats_per_plane*u.deg    #meanAnomaly steps
height = 550*u.km + const.R_earth   
incl = 53*u.deg

Pos = np.zeros([N_sats,3,steps])
indSat = 0   
fig = plt.figure(figsize=[20,10])
ax = plt.axes(projection=ccrs.PlateCarree())
ax.stock_img()
c = list()
for i in range(N_planes):
    raa = raa_step*i
    for j in range(Sats_per_plane):
        mA = mA_step*j  #mean anomaly of the satellite (i,j)
        phi = np.linspace(mA.to(u.rad).value,mA.to(u.rad).value+phi_max,steps)*u.rad # angle as a parameter
        theta = 0
         
        x = height*np.cos(phi)*np.cos(theta)
        y = height*np.cos(phi)*np.sin(theta)
        z = height*np.sin(phi)


        #rotation about X axis (equivalent to inclination)
        alpha = -(90*u.deg-incl)
        C = np.cos(alpha)
        S = np.sin(alpha)
        INC_T = np.array([[1, 0, 0],[0, C, -S],[0, S, C]])
        x1,y1,z1 = np.matmul(INC_T,np.array([x,y,z]))

        #rotation about Z axis (equivalent to ra)
        C = np.cos(raa)
        S = np.sin(raa)
        RA_T = np.array([[C, -S, 0],[S, C, 0],[0, 0, 1]])
        Pos[indSat] = np.matmul(RA_T,np.array([x1,y1,z1]))
        
        #conversion to spherical coordinaes to plot in the world map view.
        c.append(SkyCoord(Pos[indSat,0],Pos[indSat,1],Pos[indSat,2],unit='km', representation_type='cartesian', frame='gcrs', obstime = epoch))
#        ax.plot(c[-1].spherical.lon.wrap_at(180*u.deg),c[-1].spherical.lat,'o')
#        ax.plot(c[-1].spherical.lon[0].wrap_at(180*u.deg),c[-1].spherical.lat[0],'X',markersize=10)
        indSat += 1
        print (indSat)

#%% Convert the GCRS coordinates to Horzontal (or AltAz) coordinates
        
SKA_mid_core = EarthLocation.from_geodetic(lon=21.44*u.deg,lat=-30.7*u.deg,height=1000*u.m)
SKA_mid_frame = AltAz(obstime=epoch,location=SKA_mid_core)

c_AltAz = [None]*N_sats
#indVis = list()
for i in range(N_sats):
    c_AltAz[i] = c[i].transform_to(SKA_mid_frame)
#    if sum(np.array(c_AltAz[-1].alt > 0).astype('int')):
#        indVis.append(indSat)
    print (i)
    
#for i in range(N_sats):
#    np.savez('\Sat_orbits\Sat_AltAzcoords_%dsteps_%d'%(steps,i),c_AltAz = c_AltAz[i])
#    print (i)
#%% convert the list to a numpy array to operate easily
sat_pos = np.zeros([N_sats,steps,6]) #[:,:,0] = az #[:,:,1] = el #[:,:,2] = distance ; x,y,z
for i in range(N_sats):
    sat_pos[i,:,0] = c_AltAz[i].az
    sat_pos[i,:,1] = c_AltAz[i].alt
    sat_pos[i,:,2] = c_AltAz[i].distance
    sat_pos[i,:,3] = c_AltAz[i].cartesian.x
    sat_pos[i,:,4] = c_AltAz[i].cartesian.y
    sat_pos[i,:,5] = c_AltAz[i].cartesian.z
    
    print(i)
np.savez('Sat_orbits\sat_pos_%dsteps'%(steps),sat_pos = sat_pos)

#%% Plot the orbits in the AltAz frame


plt.figure()
for i in range(N_sats):
#    i = indVis[100] # uncomment to plot only one satellite
#    plt.plot(c_AltAz[i].az[c_AltAz[i].alt>=0],c_AltAz[i].alt[c_AltAz[i].alt>=0],'o') 
#    plt.plot(sat_pos[i,:,0],sat_pos[i,:,1],'o') #plot all the positions of the satellites in the AltAz frame
    plt.plot(sat_pos[i,sat_pos[i,:,1]>=0,0],sat_pos[i,sat_pos[i,:,1]>=0,1],'o')# uncomment to plot only the times that the elevation is grater than 0, (visible sat)
# Select a Alt Az to point
AltPoint = 8.2088*u.deg
AzPoint = 242.28*u.deg
plt.plot(AzPoint,AltPoint,'xr',markersize=5)



#%% plot the visible sates in a determined time

plt.figure('fig altaz')
plt.figure('fig distance')
indT = 40 #time point
count = 0
for i in range(N_sats):
    az = sat_pos[i,indT,0]
    alt = sat_pos[i,indT,1]
    if alt >0:
        count += 1
        plt.figure('fig altaz')
        plt.plot(az-180,alt,'o')
        plt.figure('fig distance')
        plt.plot(sat_pos[i,indT,2]/1e3,'o')
print(count)


#%% Calculating received power from satellites

# Creating the elev az grid
#Pointing parameters
elev_min = 15*u.deg
elev_max = 90*u.deg
elev_step = 1*u.deg
Az_min = 15*u.deg
Az_max = 90*u.deg


elev = np.linspace(elev_min.value,elev_max.value,int((elev_max.value-elev_min.value)/elev_step.value))*u.deg
N_elev = len(elev)
el = list()
az = list()

for i in range(len(elev)):
    width = np.array(90/N_elev/np.cos(elev[i]*np.pi/180)).astype('int32')
    Az = (np.linspace(-180,180,int(360/width)))*u.deg
    for j in range(len(Az)):
        el.append(elev[i].value)
        az.append(Az[j].value)

el = np.array(el)
az = np.array(az)
#%% received power
# satellite parameters
EIRP_4kHz = -15 #dBW/4kHz
BW_ch = 250*u.MHz
fmin = 10.7*u.GHz #minimum assigned frequency
fmax = 12.7*u.GHz #max assigned frequency


#calcs
freq = 11*u.GHz # np.linspace(fmin,fmax,20)
EIRP_Hz = EIRP_4kHz + 10*np.log10(1/(4*1e3))
EIRP_250MHz = EIRP_Hz + 10*np.log10(250*1e6)
EIRP_2000MHz = EIRP_Hz + 10*np.log10(2000*1e6)

EIRP_250MHz_lin = 10**(EIRP_250MHz/10)

indVis = np.where(sat_pos[:,:,1]>=0)
avePrx = np.zeros(len(el))
maxPrx = np.zeros(len(el))
Prx = np.zeros([N_sats,steps])
Prx_time = np.zeros([len(el),steps])

#el = np.array([0.705564, 1])
#az = np.array([177.663, 170])

for i in range(len(el)):
    P = np.array([np.cos(az[i]*u.deg)*np.cos(el[i]*u.deg),np.sin(az[i]*u.deg)*np.cos(el[i]*u.deg), np.sin(el[i]*u.deg)]) # pointing in cartesian coords
    blk,P = np.meshgrid(np.ones(len(indVis[0])),P)
    d = sat_pos[indVis[0],indVis[1],2]      # Distance in metres

    P_sats = np.array([sat_pos[indVis[0],indVis[1],3], sat_pos[indVis[0],indVis[1],4], sat_pos[indVis[0],indVis[1],5]])/d # satellite po in cartesian    
    
    eff_alt = np.arccos(np.einsum('ij,ij->j', P, P_sats)) * 180 / np.pi
    
    G_lin = antenna_gain_times(eff_alt)  
    FSPL_lin = (d**2) * 212521 #in times
    
    Prx[indVis[0],indVis[1]] = EIRP_250MHz_lin * G_lin/FSPL_lin
    Prx_time[i,:] = np.sum(Prx,0)
    avePrx[i] = np.sum(Prx_time[i,:])/steps
    maxPrx[i] = np.max(Prx_time[i,:])
    print(i)


#%% plot the grid on the sky
        
y = np.array(el)
x = np.array(az)
z = 10*np.log10(np.array(avePrx))
 

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

#%% plot the time domain signal in a particular direction
AltPoint = 0.7
AzPoint = 177
ind = np.where((el>=AltPoint) & (az>= AzPoint))[0][0]

plt.figure()
plt.plot(t,10*np.log10(Prx_time[ind]))
plt.title('time domain received power in Az:%f , El: %f' %(AzPoint,AltPoint))
