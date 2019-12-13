# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 22:15:32 2019

    Assesment of satellite constellation on SKA

@author: f.divruno
"""

#Commands to run in Ubuntu server:

# tmux to use a deaachable console.
# python3 to load the interactive python (so that the data remains in memory)
# subprocess.run(["git","pull"]) #to use github pull
# exec(open("./Satellite_received_power_analysis.py").read()) #to execute the script from the interactive python.





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
import os as os



# Clear the screen
os.system('clear')

if os.name == 'posix':
    print('Using Ubuntu settings')
    ubuntu = 1
if os.name == 'nt':
    print('Using Windows settings')
    ubuntu = 0
    
#ubuntu = int(input('Select Ubuntu (1) or Windows (0) : '))

if ubuntu:
# in ubuntu
    matplotlib.use('Agg')
    plt.close('all')
    outdir = r'/mnt/data/satellite_rx_power/output/'
    indir = r'/mnt/data/satellite_rx_power/input/'
else:
# in windows
    outdir = 'C:\\Users\\F.Divruno\\Dropbox (SKA)\\13- Spectrum Management\\Oneweb - SpaceX starlink\\simulation\\output\\'
    indir = 'C:\\Users\\F.Divruno\\Dropbox (SKA)\\13- Spectrum Management\\Oneweb - SpaceX starlink\\simulation\\input\\'

#%%


font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 10}

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
#    axis = np.asarray(axis)
#    axis = axis / math.sqrt(np.dot(axis, axis))
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
    epoch += np.random.random(1)*3600*u.s
    
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
        ind1 = np.where(angle<0.1)
        G[ind1] = 60
        ind2 = np.where((angle>=0.1) & (angle<=3.5))
        G[ind2] = 32-25*np.log10(angle[ind2])
        ind3 = np.where(angle>3.5)
        G[ind3] = 0
        
    else:
        if angle<0.1:
            G = 60
    
        elif (angle>0.1) & (angle<3.5):
            G = 32-25*np.log10(angle) # dB
               
        else:
            G=0
 
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
    
    Returns: 
        PSDrx_t : time vector with received total Power spectral density (i.e from all the satellites)
        
    '''
    
    N_sats = np.size(A,0)
    N_steps = np.size(A,1)
    freq = 11*u.GHz
    
    
    # transform the pointing direction to the ECEF frame
    pointingENU = np.array(Pointing_to_ENU(Pointing_el,Pointing_az))
    pointingECEF = np.matmul(T,pointingENU) #unit vector
    pointingECEF = np.reshape(np.array(pointingECEF),[3,1])

    PSDrx = np.ones([N_sats,N_steps])*-500     
   
    D = A - Rx_vect
    D_norm = np.linalg.norm(D,axis=2)
    index = np.where(D_norm<2700)
    D_norm = D_norm[index]
    [blk,Point_vect] = np.meshgrid(np.ones(len(index[0])),pointingECEF,indexing='ij')
    FSPL = Free_space((D_norm*u.km).to('m'),freq.to('MHz'))
    cos_alpha1 = ((Point_vect*D[index]).sum(1)/D_norm)
    alpha = np.arccos(cos_alpha1)*180/np.pi
    [blk, G_SKA] = antenna_gain(alpha)
    PSDrx[index] = EIRP_Hz + 30 - FSPL + G_SKA + G_sat_sidelobes #30 is to pass EIRP to dBm/Hz 
    
    PSDrx_t = np.sum(10**(PSDrx/10),0)# sum of all the satellites contributions for each timestep
        
    return PSDrx_t  

def propagate_orbits2(step,N_steps,sats,rand_seed=0,load_saved=0, trial_nr = 0):
    '''
        step: time step
        N_steps: number of time steps
        sats: list of satellite orbits
        
    '''
    if load_saved == 1:
        Aux = np.load(outdir+'Sats_propagated_'+str(step*N_steps)+'s_'+str(N_steps)+'steps_'+str(trial_nr)+'.npz')
        A = Aux['A'] 
        
    else:
            
        
        N_steps_propagate = 10
        step_propagate = N_steps*step/N_steps_propagate # coarse time step to propagate the trajectory.
        
        N_sats = len(sats)
        angle_vel = 2*np.pi/24/3600 #rad/sec
        axis = np.array([0,0,1]) # Earth obliquity is 23deg, should this be considered?    
    
        np.random.seed(int(rand_seed))
        rand_start = np.random.rand(1)*10*24*3600*u.s#time in seconds for a random start, within 10 days of T0
        theta0 = (rand_start.value)*angle_vel    # degreees of rotation of the earth
        
        A = np.zeros([N_sats,N_steps,3])    #array to contain the propagated satellites
    
        t0 = np.linspace(0,N_steps*step,N_steps_propagate)  # time array to propagate the satellites
        t1 = np.linspace(0,N_steps*step,N_steps)    #time array to interpolate the orbits with finer steps
        
        i=0    
        for sat in sats: #for each satellite
            x1 = np.zeros(N_steps_propagate)
            y1 = np.zeros(N_steps_propagate)
            z1 = np.zeros(N_steps_propagate)
            for j in range(N_steps_propagate): #propagation loop
                theta = j*step_propagate*angle_vel + theta0
                ss_delta = sat.propagate(j*step_propagate*u.s + rand_start) #propagate the orbit for the step time
                [x,y,z] = ss_delta.r # get the position in ECEF frame
                [x1[j],y1[j],z1[j]] = np.dot(rotation_matrix(axis, theta), [x.value,y.value,z.value]) # Rotate the frame along Z axis (simulate earth rotation)
                
            
            #finish propagation of the satellite, interpoate the track
            x_int = np.interp(t1,t0,x1)
            y_int = np.interp(t1,t0,y1)
            z_int = np.interp(t1,t0,z1)
            A[i,:,:] = np.transpose(np.array([x_int,y_int,z_int]))        
            i+=1
            print('Sat number %d ' %(i)) 
        np.savez_compressed(outdir + 'Sats_propagated_'+str(step*N_steps)+'s_'+str(N_steps)+'steps_'+str(trial_nr) , A=A)
    return A

#%% satellite parameters

N_planes = 24
Sats_per_plane = 66
N_sats = int(N_planes*Sats_per_plane)
raan_step = 360/N_planes        #Right ascencion steps
nu_step = 360/Sats_per_plane    #Declination steps
height = 550    
inclination = 53    



EIRP_4kHz = -15 #dBW/4kHz
BW_ch = 250*MHz
fmin = 10.7*GHz #minimum assigned frequency
fmax = 12.7*GHz #max assigned frequency
EIRP_Hz = EIRP_4kHz + 10*np.log10(1/(4*kHz))
EIRP_250MHz = EIRP_Hz + 10*np.log10(250*MHz)


G_sat_sidelobes = -10 # Gain of a sidelobe


# Receiver parameters
SKA_core_lat = -30.71278*u.deg
SKA_core_lon = 21.44305*u.deg
SKA_core_h = 1000*u.m
SKA_core_ECEF = np.array(pyproj.transform(lla,ecef,SKA_core_lon.value,SKA_core_lat.value,SKA_core_h.value))/km

SKA_core_vect = np.dot(np.reshape(np.array(SKA_core_ECEF),[3,1]),np.ones([1,N_sats]))
#SKA_core_vect = np.dot(np.reshape(np.array(SKA_core_ECEF),[3,1]),np.ones([1,N_steps]))
T = ENU_to_ECEF(SKA_core_lat,SKA_core_lon) # Transformation matrix from receiver ENU frame to ECEF frame

P_saturation = -90 #dBm


#%% Simulation control: 
# Centre frequency 
freq = 11*u.GHz # np.linspace(fmin,fmax,20)

#Number of trials
trials = 1

#Pointing parameters
elev_step = 10
elev_min = 15*u.deg
elev_max = 90*u.deg

#Number of satellites to simulate
#sat_range = np.arange(0,N_sats).astype('int') #all satellites
#sat_range = np.arange(440,442).astype('int') #in case of wanting to simulate less satellites
sat_range = [435,436,437,438,440,442] #providing a list with the number of sat


#Time control
totalTime = 6000
totalHours = totalTime/3600 #total duration of the simulation
N_steps = 6000
step = totalTime/N_steps  # time step for the propagation of the satellites.

#%% Generate the satellite constelation

sats = list()
for i in range(N_planes):
    for k in range(Sats_per_plane):
        sats.append(Generate_NGSO(height,inclination,i*raan_step,k*nu_step))
        print(i,k)
#%%
# Plot the constellation:
plot_orbits = 1
if plot_orbits:
    fig = plt.figure(figsize=[15,10])
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.stock_img()        
    _lon = np.zeros(len(sats))
    _lat = np.zeros(len(sats))
    _h = np.zeros(len(sats))
    _lon2 = np.zeros(len(sats))
    _lat2 = np.zeros(len(sats))

    for i in range(len(sats)):
        sat = sats[i]
        [_lon[i],_lat[i],_h[i]] = pyproj.transform(ecef,lla,sat.r[0].value*1e3,sat.r[1].value*1e3,sat.r[2].value*1e3)
        
    ax.plot(_lon,_lat,'-o', linewidth=1,transform=ccrs.Geodetic())
    
    print('Saving figure of satellites orbits')
    plt.xlabel('lat')
    plt.ylabel('lon')
    plt.grid()
#    plt.savefig(outdir+ 'satellites_orbits_', dpi=600, bbox_inches='tight')                        
    
#%% Calculate the grid points
elev = np.linspace(elev_min.value,elev_max.value,int((elev_max.value-elev_min.value)/elev_step))*u.deg        
N_elev = len(elev)
el = list()
az = list()

for i in range(len(elev)):
    width = np.array(90/N_elev/np.cos(elev[i])).astype('int32')
    Az = (np.linspace(-180,180,int(360/width)))*u.deg
    for j in range(len(Az)):
        el.append(elev[i])
        az.append(Az[j])

# Plot the grid:
#        TODO:

# to have discrete grid points
#el = list([15*u.deg,50*u.deg,70*u.deg])
#az= list([-105*u.deg,-105*u.deg,-105*u.deg])


#%% For each trial propagate the constellation 

k=0    
total = trials #total number of itrations
N_sats = len(sat_range) # number of satellites to simulate


A0 = np.zeros([trials,N_sats,N_steps,3]) # variable to hold the coordinates of the sallites for each trial for each timestep
for ind_trial in range(trials):
    np.random.seed(ind_trial)
    rand_seed = np.random.rand(1)*100 #generate random seed for each trial, this makes the propagation to be in a random time.
    sats2 = [sats[i] for i in sat_range]
    A0[ind_trial,:,:,:] = propagate_orbits2(step,N_steps,sats2,rand_seed,load_saved=0,trial_nr = ind_trial)  # propagate the orbits with a random seed  A= [trial, sat, time, xyz]
    print('Propagating trial: %d of %d'%(ind_trial,trials))

#A = A0[:,sat_range,:,:] # get only the number of satellites to simulate
A = A0.copy() # A is in ECEF frame, in km

#%% thansform the (XYZ) coordinates of each satellite to az elev

#Calculate the transformation matrix
T = ECEF_to_ENU(SKA_core_lat,SKA_core_lon)
sat_azeld = np.zeros([trials,N_sats,N_steps,3])
for k in range(trials):
    for i in range(N_sats):
        print(i)
        for j in range(N_steps):
            Sat_ENU = np.matmul(T,A[k,i,j,:])*u.km
            az,el = ENU_to_Pointing(Sat_ENU[0].to(u.m),Sat_ENU[1].to(u.m),Sat_ENU[2].to(u.m))            
            d = np.linalg.norm(Sat_ENU)
            sat_azeld[k,i,j,:] = [az.value,el.value,d]
            

sat_azeld[0,0,0,:]
d0 = np.linalg.norm(SKA_core_ECEF - A[0,0,0,:])

#%% Receive power 

#prepare the loop:
PSDrx = np.ones([trials,len(el),N_steps])*-500 # received power for each trial in each pointing in each time step
k=0
total = trials*len(el)     

for ind_trial in range(trials):
    [blk,blk,Rx_vect] = np.meshgrid(np.ones(N_sats),np.ones(N_steps),SKA_core_ECEF,indexing='ij')
    for ind_el in range(len(el)):  # for each pointing direction (el,az)
        PSDrx[ind_trial,ind_el,:] = received_power(Rx_vect,T,el[ind_el],az[ind_el],A[ind_trial,:,:,:]) #returns the received power for the jth pointing and ith trial for all the time steps.
        k+=1
        print('Received power completed: %.2f'%(k*100/total))    



#%% Plot orbits in the full time calculated
   # plot over the map with PLateeCarree projection.

plot_Map = 0
if plot_Map:    
    fig = plt.figure(figsize=[15,10])
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.stock_img()
    total = N_sats*trials
    k=0
    col = list(['b','r'])
    for ind_trial in range(trials): #number of trials
        for ind_sat in range(N_sats): #number of satellites to plot

            [_lon0,_lat0,_h0] = pyproj.transform(ecef,lla,A[ind_trial,ind_sat,0,0]*1e3,A[ind_trial,ind_sat,0,1]*1e3,A[ind_trial,ind_sat,0,2]*1e3)

            [_lon,_lat,_h] = pyproj.transform(ecef,lla,A[ind_trial,ind_sat,:,0]*1e3,A[ind_trial,ind_sat,:,1]*1e3,A[ind_trial,ind_sat,:,2]*1e3) #plot all the timesteps for a the ith satellite.

            ax.plot(_lon,_lat,'-o',color=col[ind_trial], linewidth=0.5,transform=ccrs.Geodetic())
            ax.scatter(_lon0,_lat0,marker='x',color='r', linewidth=5,transform=ccrs.Geodetic())
            print('Drawing trajectory, completed: %.2f'%(k/total*100))
            k+=1



#%% Calculate the average and the maximum values.
PSDrx_linear = PSDrx

PSD_ave = 10*np.log10(np.sum(np.sum(PSDrx_linear,2),0)/total) # dimmension 0 is the numbers of trials, dimmension 2 is the time steps
PSD_max = 10*np.log10(np.max(np.max(PSDrx_linear,2),0))



#%% plot the result

#z = np.array(PSD_ave+10*np.log10(BW_ch))
z = np.array(PSD_max+10*np.log10(BW_ch))

#generate the grid to plot
el2 = np.zeros(len(el))
az2 = np.zeros(len(az))
for i in range(len(el)):
    el2[i] = el[i].value        
    
for i in range(len(az)):
    az2[i] = az[i].value        
    
y = np.array(el2)
x = np.array(az2)


import matplotlib.tri as mtri

triang = mtri.Triangulation(x, y)

fig = plt.figure()
ax = fig.add_subplot(1,1,1,projection='3d')

axTris = ax.plot_trisurf(triang, z, cmap='jet')

ax.view_init(elev=90, azim=-90)

cbarAx = fig.add_axes((0.75, 0.1, 0.02, 0.8))
cbar = fig.colorbar(axTris,cax=cbarAx,orientation='vertical')

ax.set_xlabel('Azimuth')
ax.set_ylabel('Elevation')
ax.set_zlabel('Ave power dBm ')

plt.title('Received max power in %d trials, %d steps of %.2f seg'%(trials,N_steps,step))
plt.show()
print('Saving figure of received power')
plt.savefig(outdir+ 'received_power_max', dpi=600, bbox_inches='tight')                        



#%% plot the time vector in one pointing
elevation = 40
azimuth = -106
b = np.where((x>azimuth) & (y>elevation))[0]

PSD_time = PSDrx[0,b[0],:]
Pow_time = 10*np.log10(PSD_time) + 10*np.log10(BW_ch)
 
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(Pow_time)
plt.title('Received power in time - El: %d  Az: %d'%(elevation, azimuth))
plt.savefig(outdir+ 'received_power_time', dpi=600, bbox_inches='tight')  



#%% Calculate the agregated power as a function of time
# Consier case of SKA-MID pointing to zenith and  the NGSO not pointing to SKA antenna
plot_direction =0
if plot_direction:
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
    
    PSDrx = np.ones([N_steps,int(N_sats)])*-500 
    alphaList = list()
    SKA_core_norm = np.linalg.norm(SKA_core_ECEF)
    Sat_norm = np.linalg.norm(A[0,0,:]) #all satellites are the same distance from the centre of the earth.
    SKA_core_vect = np.dot(np.reshape(np.array(SKA_core_ECEF),[3,1]),np.ones([1,N_sats]))
    
    # ==================== DEBUG
    #N_steps = 1  # force the time steps to one.
    
    # ==================== DEBUG
        
    g_list1 = list()
    for i in range(N_steps):
        D = np.transpose(A[0,:,i,:])-SKA_core_vect
        D_norm = np.linalg.norm(D,axis=0)*u.km
        index = np.where(D_norm<2700*u.km)[0].astype('int32')
        
        FSPL = Free_space(D_norm[index].to('m'),freq.to('MHz'))
        
        Point_vect = np.dot(np.reshape(np.array(pointingECEF),[3,1]),np.ones([1,len(index)]))
        cos_alpha1 = ((Point_vect*D[:,index]).sum(0)/D_norm[index].value)    
    
        alpha = np.arccos(cos_alpha1)*180/np.pi
        
        alphaList.append(alpha)
        [blk, G_SKA] = antenna_gain(alpha)
    
        satLon1,satLat1,alt = pyproj.transform(ecef,lla,A[0,index,i,0]*1e3,A[0,index,i,1]*1e3,A[0,index,i,2]*1e3)
        ax.scatter(satLon1,satLat1,s=10,c = 10**(G_SKA/10),transform=ccrs.Geodetic(),cmap='jet',alpha=1)    
    #    ax.scatter(satLon1,satLat1,s=10,c =(1-alpha/90),transform=ccrs.Geodetic(),cmap='jet',alpha=1)    
    
        g_list1.append(G_SKA)
        print('Calculating received power, completed: %.2f'%(i/N_steps*100))
        PSDrx[i,index] = EIRP_Hz + 30 - FSPL + G_SKA + G_sat_sidelobes #30 is to pass EIRP to dBm/Hz 
    
    plt.title('Angle distribution- pointing: elev= %.2f  Az= %.2f '%(pointing[0].value,pointing[1].value))
    
    PSDrx_t = 10*np.log10(np.sum(10**(PSDrx/10),1)) # sum of all the satellites contributions for each timestep
        
    
    #SKA_std = RFI.SKA_EMIEMC_std(freq,freq,20,0)
    time = np.linspace(0,step*N_steps,N_steps)
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
    
