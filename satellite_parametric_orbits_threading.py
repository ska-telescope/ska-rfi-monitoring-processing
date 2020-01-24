# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 01:26:31 2019

@author: f.divruno
"""


import numpy as np
import pycraf as pycraf
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import cartopy.crs as ccrs
import astropy.constants as const
import os as os
import astropy.units as u
from numba import jit
import concurrent
import time


#%% Functions 

def ras_pattern(
        phi, diameter, wavelength, eta_a=1
        ):
    '''
    Antenna gain as a function of angular distance after `ITU-R Rec RA.1631
    <https://www.itu.int/rec/R-REC-RA.1631-0-200305-I/en>`_.

    Parameters
    ----------
    phi :  Angular distance from looking direction [deg]
    diameter :  Antenna diameter [m]
    wavelength :  Observing wavelength [m]
    eta_a :  Antenna efficiency (default: 100%)

    Returns
    -------
    gain : Antenna gain [linear]

    Notes
    -----
    - See `ITU-R Rec. RA.1631-0
      <https://www.itu.int/rec/R-REC-RA.1631-0-200305-I/en>`_ for further
      explanations and applicability of this model.
    '''
    phi = np.abs(phi)
    # eta_a = eta_a / 100.

    # the following are independent on phi, no need to compute after broadcast
    d_wlen = diameter / wavelength

    # Note, we use the version that accounts for antenna efficiency
    # see ITU radio astronomy handbook, page 50
    # gmax = 20 * np.log10(np.pi * d_wlen)
    gmax = 10 * np.log10(eta_a * (np.pi * d_wlen) ** 2)

    g1 = -1. + 15. * np.log10(d_wlen)
    phi_m = 20. / d_wlen * np.sqrt(gmax - g1)
    phi_r = 15.85 * d_wlen ** -0.6

    # note automatic broadcasting should be possible
    # _tmp = np.broadcast(phi, _diam, _wlen)
    (
        phi, d_wlen, gmax, g1, phi_m, phi_r,
        ) = np.broadcast_arrays(
        phi, d_wlen, gmax, g1, phi_m, phi_r,
        )
    gain = np.zeros(phi.shape, np.float64)

    # case 1:
    mask = (0 <= phi) & (phi < phi_m)
    gain[mask] = gmax[mask] - 2.5e-3 * (d_wlen[mask] * phi[mask]) ** 2

    # case 2:
    mask = (phi_m <= phi) & (phi < phi_r)
    gain[mask] = g1[mask]

    # case 3:
    mask = (phi_r <= phi) & (phi < 10.)
    gain[mask] = 29. - 25. * np.log10(phi[mask])

    # case 4:
    mask = (10. <= phi) & (phi < 34.1)
    gain[mask] = 34. - 30. * np.log10(phi[mask])

    # case 5:
    mask = (34.1 <= phi) & (phi < 80.)
    gain[mask] = -12.

    # case 6:
    mask = (80. <= phi) & (phi < 120.)
    gain[mask] = -7.

    # case 7:
    mask = (120. <= phi) & (phi <= 180.)
    gain[mask] = -12.

    return 10**(gain/10)


def create_orbits(N_planes=24,
                  Sats_per_plane=66,
                  orbit_incl=53*u.deg,
                  orbit_height=550*u.km,
                  orbit_period=96*u.min, 
                  total_time=1*u.h,
                  time_steps=1000, 
                  Obs_Lat=-30.7*u.deg, 
                  Obs_Lon=21.44*u.deg, 
                  Obs_height=1000*u.m,
                  plot_flag=0,
                  rand_seed=0,
                  explicit=False,
                  plane_index=[],
                  sat_index=[],
                  ):
    '''
        Creats the orbits of all the satellites and propagates them for the amount of time indicated
        then it makes an interpolation of the XYZ position of the satellite in the AltAz ref frame
        to refine the propagation without so much computational load.
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
            explicit: flag to plot different stages of the orbit generation process
            plane_index: to generate a specific plane
            sat_index: to generate a specific satellite
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
    delta = np.random.random(1)*60*60*24*u.s # RAndom number within 1 day
    epoch0 = Time("2015-01-01 00:00:00", scale="utc")+delta
    epoch = epoch0 + t
    
    phi_max = 2*np.pi*N_orbits
    
    #Orbital parameters
    N_sats = int(N_planes*Sats_per_plane)
    raa_step = 360/N_planes*u.deg        #Right ascencion steps
    mA_step = 360/Sats_per_plane*u.deg    #meanAnomaly steps
    height = orbit_height + const.R_earth   
    incl = orbit_incl
    
    if plot_flag:    
        # Plot the trajectory of the satellites  
        plt.figure('Sat_Tracks',figsize=[20,10])
        ax_sat_track = plt.axes(projection=ccrs.PlateCarree())
        ax_sat_track.stock_img()
        # Plot the trajectory of the satellites  
    
    Pos = np.zeros([N_sats,3,steps]) #cartesian position in inertial reference system
    indSat = 0 
    c = list()
    for i in range(N_planes):
        raa = raa_step*(i)
        for j in range(Sats_per_plane):
            
            # 1st step: generate an orbit with 90 deg inclination 
            # and position every satellite in the corresponding mean Anomaly
            mA = mA_step*j  #mean anomaly of the satellite number (i,j)
            phi = np.linspace(mA.to(u.rad).value,mA.to(u.rad).value+phi_max,steps)*u.rad # angle as a parameter
            #theta = 0
             
            x = height*np.cos(phi) #height*np.cos(phi)*np.cos(theta) #cos(theta)=1
            y = np.zeros(len(phi)) #height*np.cos(phi)*np.sin(theta) #is zero since theta =0
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

            if (explicit)&(j==0)&(i==5):
                fig1=plt.figure('orbit_gen',figsize=[15,9])
                ax = fig1.add_subplot(1,1,1,projection='3d')                
                ax.plot(x,y,z,'.')
                ax.plot(x1,y1,z1,'.')
                ax.plot(Pos[indSat,0],Pos[indSat,1],Pos[indSat,2],'.')
                ax.view_init(45,-90)

               

            
            #generate the Skycoord element for each satellite
            c.append(SkyCoord(Pos[indSat,0],Pos[indSat,1],Pos[indSat,2],unit='km', representation_type='cartesian', frame='gcrs', obstime = epoch))
            
            indSat += 1
            print ('\nGenerationg orbit %d of %d, sat num %d of %d '%(i+1,N_planes,j+1,Sats_per_plane), end='\r', flush=True)            

            if plot_flag:
                plt.figure('Sat_Tracks')
                #conversion to spherical coordinaes to plot in the world map view.
                ax_sat_track.plot(c[-1].spherical.lon.wrap_at(180*u.deg),c[-1].spherical.lat,'o')
 
   
    # Observatory location in geodetic coordinates:
    Obs_location = EarthLocation.from_geodetic(lon=Obs_Lon,lat=Obs_Lat,height=Obs_height)
    
    # Generate the Observatory reference frame
    Obs_frame = AltAz(obstime=epoch,location=Obs_location)
    
    # Convert the GCRS coordinates to Horizontal (also called AltAz) coordinates
    c_AltAz = [None]*N_sats
    for i in range(N_sats):
        c_AltAz[i] = c[i].transform_to(Obs_frame)
        print ('Transformation GCRS to AzAlt, sat num: ' + str(i))
        

    # convert the list "c_AltAz" to a numpy array to operate easily
    # Inerpolate the position obtained with the orbit propagator
    # Interpolate the cartesian coordinates then recalculate the polar coords
        
    steps = int(time_steps)
    t2 = np.linspace(t[0],t[-1],steps)    
    sat_pos = np.zeros([N_sats,steps,6]) #[az,alt,distance,x,y,z]
    
    @jit
    def interpolation_sat():    
        for i in range(N_sats):
            sat_pos[i,:,3] = np.interp(t2,t,c_AltAz[i].cartesian.x)
            sat_pos[i,:,4] = np.interp(t2,t,c_AltAz[i].cartesian.y)
            sat_pos[i,:,5] = np.interp(t2,t,c_AltAz[i].cartesian.z)
            sat_pos[i,:,0] = np.arctan2(sat_pos[i,:,4],sat_pos[i,:,3])*180/np.pi*u.deg
            sat_pos[i,:,1] = np.arctan2(sat_pos[i,:,5],np.sqrt(sat_pos[i,:,3]**2+sat_pos[i,:,4]**2))*180/np.pi*u.deg
            sat_pos[i,:,2] = np.sqrt(sat_pos[i,:,3]**2+sat_pos[i,:,4]**2+sat_pos[i,:,5]**2)
            print('Interpolating sat num: ' + str(i))
            
    interpolation_sat()        
    return sat_pos, c_AltAz



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



@jit
def offbeam_angle(ska_el,ska_az,sat_el,sat_az):
    # angle between pointing vector and position of the satellite in AltAz frame
    # angles need to be in radians
    # To simplify the calculation the expresion 
    # cos(alpha) = dot(A.B)/mod(A)/mod(B) is expanded in
    # the spherical coordinates of A and B = (r.cos(az).cos(el), r.sin(az).cos(el), r.sin(el))
    
    
    CEa = np.cos(ska_el)
    CAa = np.cos(ska_az)
    SEa = np.sin(ska_el)
    SAa = np.sin(ska_az)
    CEb = np.cos(sat_el)
    CAb = np.cos(sat_az)
    SEb = np.sin(sat_el)
    SAb = np.sin(sat_az)
    
    eff_alt = np.arccos(CEa*CAa*CEb*CAb+CEa*SAa*CEb*SAb+SEa*SEb)*180/np.pi
    # eff_alt = np.arccos(np.einsum('ij,ij->j', P, P_sats)) * 180 / np.pi
    return eff_alt


def receive_power_threading(chunk):
    '''
        calculates the received power in the grid of pointings from el,az
        input:
            EIRP: in dBm
            fo: centre frequency in astropy units
            sat_pos: (sat index, time step, 6)
            el: elevation array as numpy array 
            az: az arrayas numpy array 
            chunks: number of portions to thread the vector el
            chunk: index of the portion of the vector el to calculate
            
    '''
    #Radiated power
    # instead of using the EIRP and calculating the spreading loss
    # the PFD from FCC filings is considered , this takes into account the
    # shape of the beam, ie between -65 to 65 deg (0 is nadir position from 
    # the satellite) the PFD is -146 dBW/m2
    # then considering the SKA antenna gain and the frequency th PRx is calc
    #    EIRP_dBW_Hz = -49.4 #dBW/Hz
    #    EIRP_dBm_Hz = EIRP_dBW_Hz + 30 # dBm/Hz
    #    #in 250 MHz channel
    #    EIRP = EIRP_dBm_Hz + 10*np.log10(250e6)
    #    EIRP_lin = 10**(EIRP/10) # in mW

    # because   Prx = PFD*Aeff
    #           Prx = PFD*G*lda**2/(4pi)
    #           Prx = G * (PFD*lda**2/4/pi)
    #           Prx = G * Pow_const (not considering the change in freq of different channels)
    Pow_const = 10**(Pfd_CH/10)*(3e8/fo)**2/(4*np.pi)
        
    lda = (3e8/fo) #m
    D = 14.5 # m
    
    time_steps = np.size(sat_pos,1)
    
    # slice the el vector so that each thread uses a portion of it
    if len(el) >= (chunk+1)*int(np.ceil(len(el)/chunks)):
        ind_start = (chunk)*int(np.ceil(len(el)/chunks))
        ind_end = (chunk+1)*int(np.ceil(len(el)/chunks))
    else:
        ind_start = (chunk)*int(np.ceil(len(el)/chunks))
        ind_end = len(el)    


    el2 = el[ind_start:ind_end]
    az2 = az[ind_start:ind_end]
    #prepare the loop
    avePrx = np.zeros(len(el2))
    maxPrx = np.zeros(len(el2))
    print('\nElevation vector in thread %d  %d'%(chunk,len(el2)))
    
    #looks for the time steps where the elevation is >=0
    ind = np.where(sat_pos[:,:,1]>=0)
    visible_sats = np.unique(ind[0])
    
    Prx = np.zeros([len(el2),time_steps])
    
    ind_1 = ind[1]
    ind_0 = ind[0]


   
    for i in range(len(el2)):
        el_rad = el2[i]*np.pi/180
        az_rad = az2[i]*np.pi/180
        start = time.time()

        for sat_ind in visible_sats:
            # Get the time indices where the satellite is visible
            time_ind = ind_1[ind_0 == sat_ind]
            

            # satellite versor in AzAlt frame, cartesian coords
            sat_el = (sat_pos[sat_ind,time_ind,1])*np.pi/180
            sat_az = (sat_pos[sat_ind,time_ind,0])*np.pi/180
            
            eff_alt = offbeam_angle(el_rad,az_rad,sat_el,sat_az)
            
            # Linear gain
#            GdB = pycraf.antenna.ras_pattern(eff_alt*u.deg,D,lda,do_bessel=False).value
#            G_lin = 10**(GdB/10)
            G_lin = ras_pattern(eff_alt,D, lda)

            #Power received from the satellite for each time step
            # power const is (S*lda^2/(4pi))
            # Considering S = constant as the highest value from FCC filings
            Prx_sat = G_lin  #* Pow_const this constant is taken out of the loop
            
            # Accumulate the received power from each satellite in each time step
            Prx[i,time_ind] += Prx_sat
            
        Prx[i,:] *= Pow_const
        end = time.time()
        print('\nThread loop time: %.3f sec,  total thread time: %.2f min' %(end-start,(end-start)*len(el2)/60))        
    # Average for all the time calculated
    avePrx = np.sum(Prx,1)/time_steps
    
    # Maximum power in all the time considered
    maxPrx = np.max(Prx,1)


    return Prx,avePrx, maxPrx, chunk




def plot_orbit_AltAz(sat_pos, plot_all=False, AzPoint=0, AltPoint=0):
    '''
        Plots the orbits in 2d rectangular plot
        
    '''
    N_sats = np.size(sat_pos,0)
    plt.figure(figsize=[10,8])
    if plot_all==False:
        for i in range(N_sats):
            plt.plot(sat_pos[i,sat_pos[i,:,1]>=0,0],sat_pos[i,sat_pos[i,:,1]>=0,1],'.')# uncomment to plot only the times that the elevation is grater than 0, (visible sat)
    else:
        for i in range(N_sats):
            plt.plot(sat_pos[i,:,0],sat_pos[i,:,1],'.')
        plt.xlim([-180,180])
        plt.ylim([-90,90])
        
    plt.plot(AzPoint,AltPoint,'xr',markersize=5)
    plt.grid()
    plt.xlabel('Azimuth')
    plt.ylabel('Elevation')
    


def plot_visible_sats(sat_pos,indT=40,figure=[]):
    
    ind = np.where(sat_pos[:,indT,1]>=0)[0]

    if figure != []:
        plt.figure(figure)
    else:
        plt.figure()
    plt.plot(sat_pos[ind,indT,0],sat_pos[ind,indT,1],'o')
    plt.title('Visible satellites in time index : ' + str(indT) +'- number of visible sats: '+ str(len(ind)))
    plt.xlabel('Azimuth')
    plt.ylabel('Altitude')
    
#    plt.figure()
#    plt.plot(sat_pos[ind,indT,2],'o')
#    plt.title('distance to visible satellites in time index : ' + str(indT))
#    plt.ylabel('Distance in km')
    
    print('Number of visible satellites in %d is: %d'%(indT,len(ind)))
    
    
def plot_trimesh(el,az,Z,title='', view_el = 60, view_az = -90):
    '''
    
    '''
    from mpl_toolkits.mplot3d import Axes3D
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



#%% Start rhe calculation
Results_folder = 'C:/Users/f.divruno/Documents/Satellite_results'
Results_folder = 'C:/Users/F.Divruno/Dropbox (SKA)/Python_codes/satellite_results'

if __name__ == '__main__':
    constellation = 1
    if constellation == 0:
        #Starlink Constellation
        Const_name = 'Starlink - '
        orbit_incl = 53*u.deg
        orbit_height = 550*u.km 
        orbit_period = 5760.0*u.s
        N_planes = 24 
        Sats_per_plane = 66
        max_time  = 5760*u.s 
        time_steps = 5700

    elif constellation == 1:    
        # OneWeb parameters
        Const_name = 'Oneweb - '
        orbit_incl = 87.9*u.deg
        orbit_height = 1200*u.km
        orbit_period = 6600*u.s
        N_planes = 36
        Sats_per_plane = 55
        max_time  = 6600*u.s
        time_steps = 3000

    elif constellation == 2:
        #Test Constellation
        Const_name = 'Test - '
        orbit_incl = 53*u.deg
        orbit_height = 550*u.km 
        orbit_period = 5760.0*u.s
        N_planes = 24 
        Sats_per_plane = 1
        max_time  = 5760*u.s 
        time_steps = 3000

    
    CHBW = 250e6 # Hz
    global fo
    fo = 11e9
    Pfd_4kHz = -146 # dBW/m2 pfd between -65 and 65 deg, consider worst case. 
    global Pfd_CH
    Pfd_CH = Pfd_4kHz + 10*np.log10(CHBW/4e3) + 30 #dBm/m2
        
    RS = 5 # random seed for create_orbits
    
    el_min = 20*u.deg
    el_max = 90*u.deg
    el_step = 1*u.deg

    N_threads = 20
    
    # Calculating received power from satellites
    global el,az
    el,az = generate_az_el_grid(el_min, el_max, el_step)

    print('number of points in the sky: %d'%len(el))
    
    fig = plt.figure(figsize=[20,15])
    plt.plot(az,el,'.')
    plt.xlabel('Azimuth')
    plt.ylabel('Elevation')
    plt.grid()
    plt.savefig(Results_folder + '/sky grid in rectuangular plot - el step %d.png'%int(el_step.value))

    fig = plt.figure(figsize=[20,15])
    plt.polar(az*np.pi/180,(90-el)*np.pi/180,'.')
    plt.savefig(Results_folder + '/sky grid in polar plot - el step %d.png'%int(el_step.value))
        
    N_trys = 1
    avePrx = np.zeros([N_trys,len(el)])
    maxPrx = np.zeros([N_trys,len(el)])
    Prx = np.zeros([N_trys,len(el),time_steps])
    # Generate N threads to divide the sky

    
    
    for i in range(N_trys):
        np.random.seed()
        RS = int(np.random.random()*10000)
        identifier = '- %d planes - %d sats pp - seed %d'%(N_planes,Sats_per_plane,RS) # for plotting and saving
        
        # Generate orbits
        global sat_pos
        [sat_pos,c_AltAz] = create_orbits(N_planes = N_planes,
                                        Sats_per_plane = Sats_per_plane,
                                        orbit_incl = orbit_incl,
                                        orbit_height = orbit_height,
                                        orbit_period = orbit_period, 
                                        total_time = max_time,
                                        time_steps = time_steps, 
                                        Obs_Lat = -30.7*u.deg, 
                                        Obs_Lon = 21.44*u.deg, 
                                        Obs_height = 1000*u.m,
                                        plot_flag = 1,
                                        rand_seed = RS,
                                        explicit = False)

        
        # Plot the visible satellites in AzAlt reference frame
        plot_orbit_AltAz(sat_pos,plot_all=True)
        plt.title('All satellites in Horizontal ref frame')
        plt.savefig(Results_folder + '/'+Const_name+'All sats horizontal ref'+identifier+'.png')
        
        plot_orbit_AltAz(sat_pos,plot_all=False)
        plt.title('Only visible satellites in Horizontal ref frame')
        plt.savefig(Results_folder + '/'+Const_name+'Visible sats horizontal ref'+identifier+'.png')
       
#%%
        #Calculate the received power
        # use threading to accelerate the calculation
        global chunks
        chunks = N_threads
        with concurrent.futures.ThreadPoolExecutor() as executor:
            for Prx_thread, avePrx_thread, maxPrx_thread,chunk in executor.map(receive_power_threading, range(chunks)):
                # put results into correct output list
                if len(el) >= (chunk+1)*int(np.ceil(len(el)/chunks)):
                    ind_start = (chunk)*int(np.ceil(len(el)/chunks))
                    ind_end = (chunk+1)*int(np.ceil(len(el)/chunks))
                else:
                    ind_start = (chunk)*int(np.ceil(len(el)/chunks))
                    ind_end = len(el)
                               
                Prx[i,ind_start:ind_end,:] = Prx_thread+1e-10
                avePrx[i,ind_start:ind_end] = avePrx_thread+1e-10
                maxPrx[i,ind_start:ind_end] = maxPrx_thread+1e-10
                print('Thread chunk %d finished'%chunk)
       
        
        
        # plot the received power in the sky
        blk,ax = plot_trimesh(el,az,maxPrx[i,:],'Maximum received power '+ identifier ,0,180)
        plt.savefig(Results_folder + '/'+Const_name+'Max received power - full sky '+identifier+'-side.png')
        ax.view_init(90,-90)
        plt.draw()
        plt.savefig(Results_folder + '/'+Const_name+'Max received power - full sky '+identifier+'-front.png')
                
        blk,ax = plot_trimesh(el,az,avePrx[i,:], 'Average received power in %s'%(max_time)+ identifier ,0,180)
        plt.savefig(Results_folder + '/'+Const_name+'Avg received power - full sky '+identifier+'-side.png')
        ax.view_init(90,-90)
        plt.draw()
        plt.savefig(Results_folder + '/'+Const_name+'Avg received power - full sky '+identifier+'-front.png')
        
        #max in time domain:
        k = np.where(maxPrx[i,:]==np.max(maxPrx[i,:]))
        plot_rx_power_in_time(az[k],el[k],Prx[i])
        plt.savefig(Results_folder + '/'+Const_name+'Instantaneous received power - el %.2f Az %.2f'%(el[k],az[k])+identifier+'.png')
#%%        
        #save the max and averaged power
        savefile = 1
        if savefile :
            files = os.listdir(Results_folder + '')
            filename = Results_folder + '/'+Const_name+'Satellites '+identifier
            j=0
            filename2 = filename + ' - ' + str(j)
            while filename2 in files:
                j+=1
                filename2 =  filename + ' - ' + str(j)
            np.savez(filename2,sat_pos,el=el,az=az,Prx=Prx,maxPrx=maxPrx[i],avePrx=avePrx[i])

#%% CDF
        plt.figure(figsize=[15,10])
        cdf = plt.hist(10*np.log10(Prx.flatten()),1000,density=True,cumulative=True,histtype='step')
        plt.title('Cumulative dist funct of instantaneous received power')
        i = np.where(cdf[0]>0.98)[0][0]
        plt.plot([cdf[1][i],cdf[1][i]],[0,1])
        plt.savefig(Results_folder + '/'+Const_name+'Cumulative dist funct of Prx'+identifier+'.png')
        
        plt.figure(figsize=[15,10])
        cdf = plt.hist(10*np.log10(avePrx.flatten()),1000,density=True,cumulative=True,histtype='step')
        plt.title('Cumulative dist funct of average received power')
        i = np.where(cdf[0]>0.98)[0][0]
        plt.plot([cdf[1][i],cdf[1][i]],[0,1])
        plt.savefig(Results_folder + '/'+Const_name+'Cumulative dist funct of avePrx'+identifier+'.png')
        
        plt.figure(figsize=[15,10])
        cdf = plt.hist(10*np.log10(maxPrx.flatten()),1000,density=True,cumulative=True,histtype='step')
        plt.title('Cumulative dist funct of maximum received power')
        i = np.where(cdf[0]>0.98)[0][0]
        plt.plot([cdf[1][i],cdf[1][i]],[0,1])
        plt.savefig(Results_folder + '/'+Const_name+'Cumulative dist funct of maxPrx'+identifier+'.png')

#%% Histograms
        plt.figure(figsize=[15,10])
        plt.hist(10*np.log10(Prx.flatten()),1000,density=True)
        plt.title('Histogram of instantaneous received power')
        plt.savefig(Results_folder + '/'+Const_name+'Histogram of Prx'+identifier+'.png')
        
        plt.figure(figsize=[15,10])
        plt.hist(10*np.log10(avePrx.flatten()),1000,density=True)
        plt.title('Histogram of average received power')
        plt.savefig(Results_folder + '/'+Const_name+'Histogram of avePrx'+identifier+'.png')
        
        plt.figure(figsize=[15,10])
        plt.hist(10*np.log10(maxPrx),1000,density=True)
        plt.title('Histogram of maximum received power')
        plt.savefig(Results_folder + '/'+Const_name+'Histogram of maxPrx'+identifier+'.png')

#%% Inspect the full Prx matrix to see the distribution of the 98% power value in every az el point
        
#        plt.figure(figsize=[15,10])
#        Prx_98perc = np.zeros(len(el))
#        # calculate the 98% value for each az,el value
#        for i in range(len(el)):
#            [num,bins] = np.histogram(10*np.log10(Prx[0,i,:]),1000)
#            cdf = np.cumsum(num)/np.sum(num)
#            ind = np.where(cdf>0.98)[0][0]
#            Prx_98perc[i] = bins[ind]
#            print(i)
#            
        #plot the histogram of this 98% values
        plt.figure(figsize=[15,10])
        plt.hist(Prx_98perc,100,density=True)
        plt.title('Histogram of 98% values of Prx for each az,el point')
        
        plt.figure(figsize=[15,10])
        plt.hist(Prx_98perc,100,density=True,cumulative=True,histtype='step')
        plt.title('CDF of 98% values of Prx for each az,el point')

        