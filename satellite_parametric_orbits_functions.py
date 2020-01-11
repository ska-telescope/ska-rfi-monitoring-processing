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
import pycraf.antenna as antenna
import os as os

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
    G[(angle>0.1) & (angle<=3.5)] = 1.58e3*(angle[(angle>0.1) & (angle<=3.5)]**(-2.5))
    G[angle> 3.5] = 1

    return G # gain in the beam   

def antenna_gain_RA1631(angle,do_bessel=False):
    '''
        Using pycraf to implement the antenna pattern in RA.1631
    '''
    D = 14.5*u.m
    lda = 3e8/11e9*u.m
    GdB = antenna.ras_pattern(angle*u.deg,D,lda,do_bessel=do_bessel).value
    G = 10**(GdB/10)+1e-4
    return G


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
            time_steps: number of steps to compute the position of the sats in the AzAlt frame.
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
    
    #generate the random offset in time
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
            sat_pos: (sat index, time step, 6)
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
    
    Prx = np.zeros([len(el),time_steps])
    for i in range(len(el)):
        # pointing in cartesian coords [3]
        Po = np.array([np.cos(az[i]*u.deg)*np.cos(el[i]*u.deg),np.sin(az[i]*u.deg)*np.cos(el[i]*u.deg), np.sin(el[i]*u.deg)])    

   
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
#            G_lin = antenna_gain_times(eff_alt)  
            G_lin = antenna_gain_RA1631(eff_alt,do_bessel=False)  

            # FSPL in linear units
            # FSPL = 20*log10(f_MHz) + 20*log10(d_m) - 27.55 # in dB
            FSPL_lin = (d**2)*(fo**2)*0.0017579236139586931 # d in metres, fo in MHz
            
            #Power received from the satellite for each time step
            Prx_sat = EIRP_lin * G_lin / FSPL_lin
            
            
            # Accumulate the received power from each satellite in each time step
            Prx[i,time_ind] += Prx_sat

       
        
        # calculate the maximum power in each timestep, is for debugging
#        Prx_max_time = np.maximum(Prx_max_time, Prx)
                                  
        
        # Average for all the time calculated
        avePrx[i] = np.sum(Prx[i])/steps
        
        # Maximum power in all the time considered
        maxPrx[i] = np.max(Prx[i])
        
        print('Received power, point %d of %d' %(i,len(el)))    

#    return Prx, Prx_time, avePrx, maxPrx
    return Prx,avePrx, maxPrx



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
    
    max_time  = 3600*u.s #1*u.h
    time_steps = 4000
    N_planes = 24
    Sats_per_plane = 66
    RS = 5 # random seed for create_orbits
    
    el_min = 20*u.deg
    el_max = 90*u.deg
    el_step = 1*u.deg
    
    # Calculating received power from satellites
    el,az = generate_az_el_grid(el_min, el_max, el_step)
    
    N_trys = 1
    avePrx = np.zeros([N_trys,len(el)])
    maxPrx = np.zeros([N_trys,len(el)])
    for i in range(N_trys):
        np.random.seed()
        RS = int(np.random.random()*10000)
        identifier = '- %d planes - %d sats pp - seed %d'%(N_planes,Sats_per_plane,RS) # for plotting and saving
        
        # Generate orbits
        sat_pos,c_AltAz = create_orbits(N_planes=N_planes,
                                        Sats_per_plane=Sats_per_plane,
                                        rand_seed=RS,plot_flag=0,
                                        total_time=max_time,
                                        time_steps=time_steps)
        
        # Plot the orbits in the AltAz frame
        plot_orbit_AltAz(sat_pos)
        plt.title('Orbits'+identifier)
        plt.savefig('..\satellite_results\Orbits in az el '+identifier+'-side.jpg')        
        
        # plot the visible sats in a determined time
        #plot_visible_sats(sat_pos,indT=40)
        
        # Radiated power
        EIRP_4kHz = -15 #dBW/4kHz
        
        #in 250 MHz channel
        EIRP = EIRP_4kHz + 10*np.log10(250e6/4e3)
        
        #Calculate the received power
        Prx, avePrx[i], maxPrx[i] = receive_power2(EIRP,11e9*u.Hz,sat_pos,el,az)
        
        
        
        # plot the received power in the sky
        blk,ax = plot_trimesh(el,az,maxPrx[i],'Maximum received power '+ identifier ,0,180)
        plt.savefig('..\satellite_results\Max received power - full sky '+identifier+'-side.jpg')
        ax.view_init(90,-90)
        plt.draw()
        plt.savefig('..\satellite_results\Max received power - full sky '+identifier+'-front.jpg')
                
        blk,ax = plot_trimesh(el,az,avePrx[i], 'Average received power in %s'%(max_time)+ identifier ,0,180)
        plt.savefig('..\satellite_results\Avg received power - full sky '+identifier+'-side.jpg')
        ax.view_init(90,-90)
        plt.draw()
        plt.savefig('..\satellite_results\Avg received power - full sky '+identifier+'-front.jpg')
        
        #max in time domain:
        k = np.where(maxPrx[i]==np.max(maxPrx[i]))
        plot_rx_power_in_time(az[k],el[k],Prx)
        plt.savefig('..\satellite_results\Instantaneous received power - el %.2f Az %.2f'%(el[k],az[k])+identifier+'.jpg')
        
        #save the max and averaged power
        files = os.listdir('..\satellite_results')
        filename = '..\satellite_results\Satellites '+identifier
        j=0
        filename2 = filename + ' - ' + str(j)
        while filename2 in files:
            j+=1
            filename2 =  filename + ' - ' + str(j)
        np.savez(filename2,el=el,az=az,Prx=Prx,maxPrx=maxPrx[i],avePrx=avePrx[i])
