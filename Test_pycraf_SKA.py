#!/usr/bin/env python
# coding: utf-8


from collections import namedtuple, OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from pycraf import pathprof
from pycraf import conversions as cnv
import xlrd
import pandas as pd

pathprof.SrtmConf.set(download='missing', server='viewpano')

#%%
'''
===========================================
    Coordinates
    Set the coordinates of the transmitter.
===========================================
'''


culpritsCsvFile = r'C:\Users\F.Divruno\Dropbox (SKA)\Python_codes\SKA1_Mid_EMC_analysis_locations.csv'
aux = pd.read_csv(culpritsCsvFile,comment='#')
print('\n\n\n\n')
print(aux.location)
loc_num = int(input('\n\nSelect location to analyze: '))
location = aux.location[loc_num]

Transmitter = location
lon_tx = aux.lon[loc_num]*u.deg
lat_tx = aux.lat[loc_num]*u.deg

#SKA_MID coordinates for the KAPB
#site, lon_tx, lat_tx =  'KAPB', 21.43138889* u.deg, -30.75277778 * u.deg

#SKA_MID coordinates for the EOC
#site, lon_tx, lat_tx =  'EOC', 21.993900* u.deg, -30.973975 * u.deg

#SADT shed 4 in MID location
#lon_tx, lat_tx =  21.812403* u.deg,  -30.866074 * u.deg

#SADT shed 002 in MID location
#site, lon_tx, lat_tx =  'Mid',  21.463920* u.deg,  -30.652300 * u.deg


# SKA_LOW coordinates for the CPF at 400m
#site, lon_tx, lat_tx = 'Low', 116.747085* u.deg, -26.834128 * u.deg

# SKA_LOW coordinates for the CPF at 800m (aprox)
#site, lon_tx, lat_tx = 'Low', 116.743729* u.deg, -26.833987 * u.deg

#%%
'''
===========================================
     Get the map height information
     
===========================================
'''


map_size_lon, map_size_lat = .5 * u.deg, .5 * u.deg
map_resolution = 0.001* u.deg
hprof_step = 100 * u.m

# get the data for the map
lons, lats, heightmap = pathprof.srtm_height_map(
    lon_tx, lat_tx,
    map_size_lon, map_size_lat,
    map_resolution=map_resolution,
    hprof_step = hprof_step
    )


'''
===========================================
    Sites of particular interest
     
===========================================
'''

Site = namedtuple('site', ['name', 'coord', 'pixoff', 'color'])
sites = OrderedDict([
    # ID: tuple(Name, (lon, lat), Type, height, diameter, (xoff, yoff), color)
    #('Tx', Site('Tx', (u.Quantity(lon_tx).value, u.Quantity(lat_tx).value), (20, +30), 'k')),
    #('Low_E4', Site('Low_E4', (116.75, -26.835), (60, -20), 'k')),
    ('Mid_114', Site('Mid_114', (21.452432,-30.739759), (30, -40), 'k')),
    ('Mid_121', Site('Mid_121', (21.406843,-30.803121), (-50, +20), 'k')),
#    ('Mid_124', Site('Mid_124', (21.804492, -30.84449), (-35, +25), 'k')),
#    ('Mid_020', Site('Mid_020', (21.4473, 	-30.6631), (-35, +25), 'k')),
#    ('Mid_120', Site('Mid_120', (21.6973, 	-30.8837), (-35, +25), 'k')),
#    ('Mid_124', Site('Mid_124', (21.8045, 	-30.8445), (-35, +25), 'k')),
#    ('Mid_021', Site('Mid_021', (21.9232, 	-30.7791), (-35, +25), 'k')),
#    ('Mid_119', Site('Mid_119', (21.3971606, 	-30.7785514), (-35, +25), 'k')),        
    ])

#%%
'''
===========================================
    Antennas locations
     
===========================================
'''

Antenna_locations = xlrd.open_workbook(r'C:\Users\F.Divruno\Dropbox (SKA)\Python_codes\SKA1_Mid_coordinates.xlsx')

sheet_indx = 0
sheet = Antenna_locations.sheet_by_index(sheet_indx)
name_x = list()
lat_x = np.zeros(sheet.nrows-2)
long_x = np.zeros(sheet.nrows-2)
for i in range(2,sheet.nrows-1):
        name_x.append(sheet.cell_value(i,0))
        long_x[i-2] = sheet.cell_value(i,1)
        lat_x[i-2] = sheet.cell_value(i,2)
        
#%%
'''
===========================================
    Plot terrain information
     
===========================================
'''

_lons = lons.to(u.deg).value
_lats = lats.to(u.deg).value
_heightmap = heightmap.to(u.m).value

fig = plt.figure(figsize=(10, 10))
ax = fig.add_axes((0., 0., 1.0, 1.0))

vmin,vmax = np.min(_heightmap),np.max(_heightmap) 

terrain_cmap,terrain_norm = pathprof.terrain_cmap_factory(sealevel=vmin,vmax=vmax)

cim = ax.imshow(
                _heightmap,
                origin='lower', interpolation='nearest',
                cmap=terrain_cmap, norm=terrain_norm,
                vmin=vmin, vmax=vmax,
                extent=(_lons[0], _lons[-1], _lats[0], _lats[-1]),
                )



ax.set_xlabel('Longitude [deg]')
ax.set_ylabel('Latitude [deg]')


# Annotate the sites of interest in the map
for site_id, site in sites.items():
    color = site.color
    color = 'k'
    ax.annotate(
        site.name, xy=site.coord, xytext=site.pixoff, 
        textcoords='offset points', color=color,
        arrowprops=dict(arrowstyle="->", color=color)
        )
    ax.scatter(site.coord[0], site.coord[1],marker='o', c='k')

# set aspect ratio and limits of the image
ax.set_aspect(abs(_lons[-1] - _lons[0]) / abs(_lats[-1] - _lats[0])) 
ax.set_ylim([_lats[0],_lats[-1]])
ax.set_xlim([_lons[0],_lons[-1]])


# Place a marker in each of the SKA antennas in the map
for i in range(len(long_x)):
    ax.scatter(long_x[i], lat_x[i],marker='o', c='k')



'''
===========================================
    Calculate the attenuation map
     
===========================================
'''
# Here input the patameters for the ITU-R 452-16 model

freq_string = input('\nSpecify the frequency range (in GHz): ')
freq = int(freq_string) * u.GHz
omega = 0. * u.percent  # fraction of path over sea
temperature = 290. * u.K
pressure = 1013. * u.hPa
timepercent = 0.02 * u.percent  # see P.452 for explanation
h_tg, h_rg = 3 * u.m, 2 * u.m#height of the receiver and transmitter above gnd
G_t, G_r = 0 * cnv.dBi, 0 * cnv.dBi
zone_t, zone_r = pathprof.CLUTTER.UNKNOWN, pathprof.CLUTTER.UNKNOWN
hprof_step = 10 * u.m


## Fasster approach

hprof_cache = pathprof.height_map_data(lon_tx, lat_tx,
                                        map_size_lon, map_size_lat,
                                        map_resolution=map_resolution,
                                        zone_t=zone_t, zone_r=zone_r,
                                        )

results = pathprof.atten_map_fast(freq,
                                temperature,
                                pressure,
                                h_tg, h_rg,
                                timepercent,
                                hprof_cache,  # dict_like
                                )


_lons = hprof_cache['xcoords']
_lats = hprof_cache['ycoords']
_total_atten = results['L_b']  # L_b is the total attenuation, considering all the factors.
_fspl_atten = results['L_bfsg']  # considers only the free space loss
 


'''
===========================================
    Plot the resulted attenuation map
    
===========================================
'''

# Plot the results selected
vmin, vmax = 90, 200 # Max and min scale

fig = plt.figure(figsize=(10, 10))
ax = fig.add_axes((0., 0.1, 1.0, 0.8))
cbax = fig.add_axes((0.3, 0., 0.4, .02))

cim = ax.imshow(
                _total_atten.to(cnv.dB).value,
                origin='lower', interpolation='nearest', cmap='inferno_r',
                vmin=vmin, vmax=vmax,
                extent=(_lons[0], _lons[-1], _lats[0], _lats[-1]),
                )

cbar = fig.colorbar(
                cim, cax=cbax, orientation='horizontal',
                )

ax.set_aspect(abs(_lons[-1] - _lons[0]) / abs(_lats[-1] - _lats[0]))

cbar.set_label(r'Path propagation loss', color='w')
cbax.xaxis.set_label_position('top')
for t in cbax.xaxis.get_major_ticks():
    t.tick2On = True
    t.label2On = True
ctics = np.arange(vmin, vmax, 10)
cbar.set_ticks(ctics)

ax.set_xlabel('Longitude [deg]')
ax.set_ylabel('Latitude [deg]')
ax.set_autoscale_on(False)
cbax.set_autoscale_on(False)


# Annotate the coordinates of the interest sites and searches for the attenuation levels.

lat_mesh, lon_mesh = np.meshgrid(_lats,_lons) # hace un mesh para buscar los puntos
for site_id, site in sites.items():
    color = site.color
    color = 'b'
    aux = abs(lat_mesh-site.coord[1])+abs(lon_mesh-site.coord[0])
    i,j = np.unravel_index(aux.argmin(),aux.shape)
    ax.annotate(
        site.name + ' att: ' + str(_total_atten.to(cnv.dB).value[j,i])[0:6] + ' dB', xy=site.coord, xytext=site.pixoff, 
        textcoords='offset points', color=color,
        arrowprops=dict(arrowstyle="->", color=color)
        )
    ax.scatter(site.coord[0], site.coord[1],marker='o', c='b') 

for i in range(len(long_x)): # Puts a mark in every place there is an antenna
    ax.scatter(long_x[i], lat_x[i],marker='o', c='w')

plt.title('Attenuation map, Freq = '+ str(freq) + ' GHz')
    
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
plt.show()

#%%

'''
===========================================
    Plot the terrain map with attenuation contours
    Also add the contours for Free Space Path Loss. 
    
===========================================
'''

_lons = lons.to(u.deg).value
_lats = lats.to(u.deg).value
_heightmap = heightmap.to(u.m).value

fig = plt.figure(figsize=(10, 10))
ax = fig.add_axes((0., 0.1, 1.0, 0.8))

vmin,vmax = np.min(_heightmap),np.max(_heightmap) 

terrain_cmap,terrain_norm = pathprof.terrain_cmap_factory(sealevel=vmin,vmax=vmax)

cim = ax.imshow(
    _heightmap,
    origin='lower', interpolation='nearest',
    cmap=terrain_cmap, norm=terrain_norm,
    vmin=vmin, vmax=vmax,
    extent=(_lons[0], _lons[-1], _lats[0], _lats[-1]),
    )

_fspl_atten = results['L_bfsg'] 

ax.contour(_total_atten.to(cnv.dB).value, levels=[100,120,140,160],
           linestyles='-',
           origin='lower',
           extent=(_lons[0], _lons[-1], _lats[0], _lats[-1]),
           alpha=1)

ax.contour(_fspl_atten.to(cnv.dB).value, levels=[100],
           colors=['red'], linestyles='-',
           origin='lower',
           extent=(_lons[0], _lons[-1], _lats[0], _lats[-1]),
           alpha=1)


'''
ax.contourf(_total_atten.to(cnv.dB).value, levels=[0,90],
           colors='b', linestyles='-',
           origin='lower',
           extent=(_lons[0], _lons[-1], _lats[0], _lats[-1]),
           alpha=0.2)

ax.contourf(_total_atten.to(cnv.dB).value, levels=[120,130],
           colors='r', linestyles='-',
           origin='lower',
           extent=(_lons[0], _lons[-1], _lats[0], _lats[-1]),
           alpha=0.4)

'''
# Position the SKA antennas in the map
for i in range(len(long_x)):
    ax.scatter(long_x[i], lat_x[i],marker='o', c='k')      
    
ax.set_ylim([_lats[0],_lats[-1]])
ax.set_xlim([_lons[0],_lons[-1]])



'''
===========================================
     Find the attenuation at each position of SKA antennas usign the map info
===========================================
'''

lat_mesh, lon_mesh = np.meshgrid(_lats,_lons) # mesh in lats and longs

for k in range(len(long_x)-1):
    
    aux = abs(lat_mesh-lat_x[k])+abs(lon_mesh-long_x[k])
    i,j = np.unravel_index(aux.argmin(),aux.shape)
    print('Antenna: %s - Att: %.2f'%(name_x[k],(_total_atten.to(cnv.dB).value[j,i])))

#%%

'''
===========================================
    Find the attenuation at each position of SKA antennas using path info
===========================================
'''
if input('find attenuation at each antenna? (Y - N)')=='Y':
    N_ant = len(long_x)-1
    N_freqs = 50
    freqs = np.logspace(np.log10(350),np.log10(15400),N_freqs)*u.MHz
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
    Atten_ant = np.zeros([N_ant,N_freqs])
    results=[]
    hprof_cache =[]
    hprof_step=100*u.m
    
    for k in range(N_ant):
        hprof_cache = pathprof.height_path_data(lon_tx,
                                                lat_tx,
                                                long_x[k]*u.deg,
                                                lat_x[k]*u.deg,
                                                hprof_step)
        for i in range(N_freqs):
            results = pathprof.atten_path_fast(freqs[i],
                                            temperature,
                                            pressure,
                                            h_tg, h_rg,
                                            timepercent,
                                            hprof_cache,  # dict_like
                                            )
            Atten_ant[k,i] = results['L_b'][-1].value #gets the last value of the attenuation path.
    #        _total_atten = results['L_b']  # L_b is the total attenuation, considering all the factors.
    #        _fspl_atten = results['L_bfsg']  # considers only the free space loss
            print(i,k,' atten value ',Atten_ant[k,i]) 
            
    fig = plt.figure(figsize=(10, 10))
    plt.semilogx(freqs,np.transpose(Atten_ant))
    plt.grid()
