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


def create_orbits(N_planes=24, Sats_per_plane=66, orbit_incl=53*u.deg,
                  orbit_height=550*u.km, orbit_period=96*u.min, total_time=1*u.h,
                  time_steps=1000, Obs_Lat=-30.7*u.deg, 
                  Obs_Lon=21.44*u.deg, Obs_height=1000*u.m,
                  plot_flag=0,
                  rand_seed=0):
    '''
        Creats the orbits of all the satellites and propagates them for the amount of time indicated
        then it makes an interpolation of the XYZ position of the satellite in the AltAz ref frame
        to refine the propagation without so much computational overhead.
        Inputs: (all quantities from astropy.quantity)
            N_planes : number or orbital planes
            Sats_per_plane: num of satellites in an orbital plane
            orbit_incl: inclination of the orbit wrt the equator in u.deg
            orbit_height: height of the orbits in u.m or u.km
            orbit_period:period in u.s
            total_time: total time of the orbit generation in u.s
            time_steps: number of steps to compute the position of the sates in the AzAlt frame.
            Obs_Lat: observer latitude
            Obs_Lon: observer longitude
            Obs_height: observer altitude wrt sea level.
            rand_seed: random seed to generate the offset in time
        output:
            sat_pos: coordinates in AltAz coord system
            c_AltAz: Original Skycoord elements with the calculated positions.
    '''
    period = orbit_period 
    steps = (int(total_time/orbit_period)+1)*100 #100 steps as a standard for the first orbit generation.
    if steps == 0:
        raise SystemExit
        
    t = np.linspace(0,1,steps)*total_time
    N_orbits = (total_time.to(u.s)/period.to(u.s)).value
    
    #generate the randome offset in time
    np.random.seed(rand_seed)
    delta = np.random.random(1)*60*60*24*u.s
    epoch0 = Time("2015-01-01 00:00:00", scale="utc")+delta
    epoch = epoch0 + t
    
    phi_max = 2*np.pi*N_orbits
    
#    N_planes = 24
#    Sats_per_plane = 66
    N_sats = int(N_planes*Sats_per_plane)
    raa_step = 360/N_planes*u.deg        #Right ascencion steps
    mA_step = 360/Sats_per_plane*u.deg    #meanAnomaly steps
    height = orbit_height + const.R_earth   
    incl = orbit_incl
    
    if plot_flag:    
        # Plot the trajectory of the satellites  
        plt.figure(figsize=[20,10])
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.stock_img()
        # Plot the trajectory of the satellites  
    
    Pos = np.zeros([N_sats,3,steps])
    indSat = 0 
    c = list()
    for i in range(N_planes):
        raa = raa_step*i
        for j in range(Sats_per_plane):
            
            # position every satellite in the corresponding mean Anomaly
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
    
            #rotation about Z axis (equivalent to right ascension)
            C = np.cos(raa)
            S = np.sin(raa)
            RA_T = np.array([[C, -S, 0],[S, C, 0],[0, 0, 1]])
            Pos[indSat] = np.matmul(RA_T,np.array([x1,y1,z1]))
            
            #generate the Skycoord element for each satellite
            c.append(SkyCoord(Pos[indSat,0],Pos[indSat,1],Pos[indSat,2],unit='km', representation_type='cartesian', frame='gcrs', obstime = epoch))
            
            indSat += 1
            print ('Generationg orbit, sat num: ' + str(indSat))
            
            if plot_flag:
                #conversion to spherical coordinaes to plot in the world map view.
                ax.plot(c[-1].spherical.lon.wrap_at(180*u.deg),c[-1].spherical.lat,'o')
 
    # Convert the GCRS coordinates to Horzontal (or AltAz) coordinates
    Obs_location = EarthLocation.from_geodetic(lon=Obs_Lon,lat=Obs_Lat,height=Obs_height)
    Obs_frame = AltAz(obstime=epoch,location=Obs_location)
    c_AltAz = [None]*N_sats
    for i in range(N_sats):
        c_AltAz[i] = c[i].transform_to(Obs_frame)
        print ('Transformation GCRS to AzAlt, sat num: ' + str(i))
        

    # convert the list to a numpy array to operate easily
    # Inerpolate the position obtained with the orbit propagator
    # Interpolate the cartesian coordinates then recalculate the polar coords
        
    steps = int(time_steps)
    t2 = np.linspace(t[0],t[-1],steps)    
    sat_pos = np.zeros([N_sats,steps,6]) #[az,alt,distance,x,y,z]
    for i in range(N_sats):
        sat_pos[i,:,3] = np.interp(t2,t,c_AltAz[i].cartesian.x)
        sat_pos[i,:,4] = np.interp(t2,t,c_AltAz[i].cartesian.y)
        sat_pos[i,:,5] = np.interp(t2,t,c_AltAz[i].cartesian.z)
        sat_pos[i,:,0] = np.arctan2(sat_pos[i,:,4],sat_pos[i,:,3])*180/np.pi*u.deg
        sat_pos[i,:,1] = np.arctan2(sat_pos[i,:,5],np.sqrt(sat_pos[i,:,3]**2+sat_pos[i,:,4]**2))*180/np.pi*u.deg
        sat_pos[i,:,2] = np.sqrt(sat_pos[i,:,3]**2+sat_pos[i,:,4]**2+sat_pos[i,:,5]**2)
        print('Interpolating sat num: ' + str(i))
        
    return sat_pos, c_AltAz


def plot_orbit_AltAz(sat_pos, AzPoint=0, AltPoint=0):
    '''
        Plots the orbits in 2d rectangular plot
        
    '''
    N_sats = np.size(sat_pos,0)
    plt.figure()
    for i in range(N_sats):
        plt.plot(sat_pos[i,sat_pos[i,:,1]>=0,0],sat_pos[i,sat_pos[i,:,1]>=0,1],'o')# uncomment to plot only the times that the elevation is grater than 0, (visible sat)
    
    plt.plot(AzPoint,AltPoint,'xr',markersize=5)


def plot_visible_sats(sat_pos,indT=40,):
    
    ind = np.where(sat_pos[:,indT,1]>=0)[0]

    plt.figure()
    plt.plot(sat_pos[ind,indT,0],sat_pos[ind,indT,1],'o')
    plt.title('visible satellites in time index : ' + str(indT))
    plt.xlabel('Azimuth')
    plt.ylabel('Altitude')
    
    plt.figure()
    plt.plot(sat_pos[ind,indT,2],'o')
    plt.title('distance to visible satellites in time index : ' + str(indT))
    plt.ylabel('Distance in km')
    
    print('Number of visible satellites in %d is: %d'%(indT,len(ind)))


def generate_az_el_grid(el_min=15*u.deg,el_max=90*u.deg,el_step=1*u.deg):
    '''
        Generates the el az grid according to the ITU-R recommendation on 
        epfd calculation.
    '''
 
    # Creating the elev az grid
    #Pointing parameters
    elev_min = el_min
    elev_max = el_max
    elev_step = el_step
     
    elev = np.linspace(elev_min.value,elev_max.value,int((elev_max.value-elev_min.value)/elev_step.value))*u.deg
    N_elev = len(elev)
    el = list()
    az = list()
    
    for i in range(len(elev)):
        width = (elev[-1] - elev[0]).value/N_elev/np.cos(elev[i]*np.pi/180)
        Az = (np.linspace(-180,180,int(360/width)))*u.deg
        for j in range(len(Az)):
            el.append(elev[i].value)
            az.append(Az[j].value)
    
    el = np.array(el)
    az = np.array(az)  
    return el, az


def receive_power2(EIRP,fo,sat_pos,el,az):
    '''
        calculates the received power in the grid of pointings from el,az
        input:
            EIRP: in dBm
            fo: centre frequency in astropy units
            sat_pos:
            el: elevation array as numpy array 
            az: az arrayas numpy array 
            
    '''
    EIRP_lin = 10**(EIRP/10)
    fo = fo.to(u.MHz).value
    
    steps = np.size(sat_pos,1)
    
    #prepare the loop
    avePrx = np.zeros(len(el))
    maxPrx = np.zeros(len(el))

    
    #To test with only one pointing
    #el = np.array([ 32.7832])
    #az = np.array([ 147.152])
    
    ind = np.where(sat_pos[:,:,1]>=0)
    visible_sats = np.unique(ind[0])
    
    for i in range(len(el)):
        # pointing in cartesian coords [3]
        Po = np.array([np.cos(az[i]*u.deg)*np.cos(el[i]*u.deg),np.sin(az[i]*u.deg)*np.cos(el[i]*u.deg), np.sin(el[i]*u.deg)])    

        Prx = np.zeros(steps)    
    
        for sat_ind in visible_sats:
            # Get the time indices where the satellite is visible
            time_ind = ind[1][ind[0] == sat_ind]
            
            # meshgrid to match the shape of the visible satellites            
            blk,P = np.meshgrid(np.ones(len(time_ind)),Po)
            
            # Distance in metres
            d = sat_pos[sat_ind,time_ind,2]      
        
            # satellite versor in AzAlt frame, cartesian coords
            P_sats = np.array([sat_pos[sat_ind,time_ind,3], sat_pos[sat_ind,time_ind,4], sat_pos[sat_ind,time_ind,5]])/d    
        
         
            #angle between pointing vector and position of the satellite in AltAz frame
            eff_alt = np.arccos(np.einsum('ij,ij->j', P, P_sats)) * 180 / np.pi
            
            # Linear gain
            G_lin = antenna_gain_times(eff_alt)  
            
            # FSPL in linear units
            # FSPL = 20*log10(f_MHz) + 20*log10(d_m) - 27.55 # in dB
            FSPL_lin = (d**2)*(fo**2)*0.0017579236139586931 # d in metres, fo in MHz
            
            # Accumulate the received power from all the satellites in each time step
            Prx[time_ind] += np.sum(EIRP_lin * G_lin / FSPL_lin,0)

       
        # calculate the maximum power in each timestep, is for debugging
#        Prx_max_time = np.maximum(Prx_max_time, Prx)
                                  
        
        # Average for all the time calculated
        avePrx[i] = np.sum(Prx)/steps
        
        # Maximum power in all the time considered
        maxPrx[i] = np.max(Prx)
        
        print('Received power, point %d of %d' %(i,len(el)))    

#    return Prx, Prx_time, avePrx, maxPrx
    return avePrx, maxPrx



def plot_trimesh(el,az,Z,title=''):
    '''
    
    '''
    
    import matplotlib.tri as mtri                
    y = np.array(el)
    x = np.array(az)
    z = 10*np.log10(np.array(Z))
     
    triang = mtri.Triangulation(x, y)
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1,projection='3d')
    
    ax.plot_trisurf(triang, z, cmap='jet', edgecolor='none', linewidth=0, antialiased=False)
    
    ax.view_init(elev=60, azim=-90)
    
    ax.set_xlabel('Azimuth')
    ax.set_ylabel('Elevation')
    ax.set_zlabel('Averaged power dBm ')
    plt.title(title)
    return fig
    

def plot_rx_power_in_time(Point_az, Point_el, Prx_time,fig=[]):
    
    ind = np.where((el>=Point_el) & (az>= Point_az))[0][0]
    
    if fig==[]:
        plt.figure()
    else:
        plt.figure(fig)
    plt.plot(10*np.log10(Prx_time[ind]))
#    plt.title('time domain received power in Az:%f , El: %f' %(Point_az,Point_el))


#%% Start rhe calculation

max_time  = 1*u.h
time_steps = 2000
N_planes = 5
Sats_per_plane = 1

el_min = 10*u.deg
el_max = 90*u.deg
el_step = 5*u.deg

# Generate orbits
sat_pos,c_AltAz = create_orbits(N_planes=N_planes, Sats_per_plane=Sats_per_plane,rand_seed=1,plot_flag=1,total_time=max_time,time_steps=time_steps)

# Plot the orbits in the AltAz frame
plot_orbit_AltAz(sat_pos)


# plot the visible sates in a determined time
#plot_visible_sats(sat_pos,indT=40)

# Calculating received power from satellites
el,az = generate_az_el_grid(el_min=el_min,el_max=el_max,el_step=el_step)

# received power
EIRP_4kHz = -15 #dBW/4kHz
#in 250 MHz channel
EIRP = EIRP_4kHz + 10*np.log10(250e6/4e3)
#Prx, avePrx, maxPrx, Prx_max_time = receive_power(EIRP,11e9*u.Hz,sat_pos,el,az)
avePrx, maxPrx = receive_power2(EIRP,11e9*u.Hz,sat_pos,el,az)

# plot the grid on the sky
fig_max = plot_trimesh(el,az,maxPrx,'Maximum received power')
fig_ave = plot_trimesh(el,az,avePrx, 'Average received power in %d s'% max_time.to(u.s).value)

