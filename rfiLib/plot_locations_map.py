# -*- coding: utf-8 -*-

"""
Created on Wed Aug 14 08:57:56 2019
    Plot dots in lon,lat over a image of the earth
    
    if pycraf is installed can use high resolution data to plot the image
    
    
@author: f.divruno
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import astropy.units as u
import numpy as np
import pandas as pd
import os
    
def plot_locations_map(lonData,latData,tx_name,centre=[],skaAntPosCsvFile='',mapsize=[5,5],scatter=1,outdir=''):
    '''
        description: plots the locations of culprits on a map and the locations of SKA anennas.
        inputs:
            centre: [lon,lat] centre of the map
            lonData, latData: culprits locations, numpy arrays in astropy units u.deg
            skaAntPosCsvFile: csv file containing the columns name,lon,lat of the SKA antennas
            mapsize: size of the map in degrees in astropy units
            
    '''
    try:
        from pycraf import pathprof
        from pycraf.pathprof import SrtmConf
        
        SrtmConf.set(download='missing', server='viewpano')
        try:
            os.mkdir('SRTMdata_temp')
        except:
            blk=1
            
        SrtmConf.set(srtm_dir='SRTMdata_temp')
        
        have_pycraf = 1
        print('Pycraf installed, plotting high res map')
    except:
        print('Pycraf not installed, plotting low res')
        have_pycraf = 0
    
    plotSkaAnt = 0
    if skaAntPosCsvFile!='':
        skaAntPos = pd.read_csv(skaAntPosCsvFile,comment='#')
        plotSkaAnt = 1

    if len(centre) == 0: #if centre is not provided calculate the average from the data
        lonCent = np.average(lonData.to(u.deg).value)*u.deg
        latCent = np.average(latData.to(u.deg).value)*u.deg
    else:
        lonCent = centre[0]*u.deg 
        latCent = centre[1]*u.deg 
    
#    have_pycraf = 0
    if have_pycraf:
        #get the heightmap information        
        map_size_lon, map_size_lat = mapsize[0] * u.deg, mapsize[1] * u.deg
        map_resolution = 0.01* u.deg
        hprof_step = 10* u.m

            
        # get the data for the map
        lons, lats, heightmap = pathprof.srtm_height_map(
                                                        lonCent, latCent,
                                                        map_size_lon, map_size_lat,
                                                        map_resolution=map_resolution,
                                                        hprof_step = hprof_step
                                                        )
        
        _lons = lons.to(u.deg).value
        _lats = lats.to(u.deg).value
        _heightmap = heightmap.to(u.m).value
        
        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
        
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
        for i in range(len(tx_name)):
            color = 'k'
            ax.annotate(
                tx_name[i], xy=np.array([lonData[i].to(u.deg).value,latData[i].to(u.deg).value]), xytext=[2,2], 
                textcoords='offset points', color=color,
                )
            ax.scatter(lonData[i], latData[i],marker='o', c='r')
        
        # set aspect ratio and limits of the image
        ax.set_aspect(abs(_lons[-1] - _lons[0]) / abs(_lats[-1] - _lats[0])) 
        ax.set_ylim([_lats[0],_lats[-1]])
        ax.set_xlim([_lons[0],_lons[-1]])
        
        
        # Place a marker in each of the SKA antennas in the map
#        for i in range(len(lonData)):
#            ax.scatter(lonData[i], latData[i],marker='o', c='r')
            
        if plotSkaAnt:
            lonSka = np.array(skaAntPos['lon'])
            latSka = np.array(skaAntPos['lat'])
            ax.scatter(lonSka, latSka,marker='o', c='lightgreen')
            
        
    else:
        lons = lonData.to(u.deg).value
        lats = latData.to(u.deg).value
        
        fig = plt.figure(figsize=[15,10])
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.stock_img()

        if plotSkaAnt:
            lonSka = np.array(skaAntPos['lon'])
            latSka = np.array(skaAntPos['lat'])
            ax.scatter(lonSka, latSka,marker='o', c='green')


        # Annotate the sites of interest in the map
        for i in range(len(tx_name)):
            color = 'k'
            ax.annotate(
                tx_name[i], xy=np.array([lonData[i].to(u.deg).value,latData[i].to(u.deg).value]), xytext=[2,2], 
                textcoords='offset points', color=color,
                )
        ax.scatter(lons,lats,marker='o',color='r', linewidth=3,transform=ccrs.Geodetic())        

        plt.xlim(np.array([-0.5*mapsize[0],0.5*mapsize[0]])+lonCent.to(u.deg).value)
        plt.ylim(np.array([-0.5*mapsize[1],0.5*mapsize[1]])+latCent.to(u.deg).value)
        
        if outdir!='':
            plt.savefig(outdir+ 'map_plot', dpi=600, bbox_inches='tight')                        
        


if __name__=='__main__':
    # autotest
    # provide the location of the SKA antennas or remove the skaAntPosCsvFile parameter.
    
    lons = np.array([22.2285    , 23.0498   , 20.604202  , 21    , 20.0351    , 21.548152, 22.270415])*u.deg
    lats = np.array([-31.0737   , -30.4068  ,-30.292972 , -31   , -31.3665   , -30.845583, -30.561479])*u.deg

    skaAntPosCsvFile = r'C:\Users\F.Divruno\Dropbox (SKA)\Python_codes\SKA1_Mid_coordinates.csv'
    
    centre = [np.average(lons.to(u.deg).value)*u.deg, np.average(lats.to(u.deg).value)*u.deg]
    centre = [21.45*u.deg,-30.7*u.deg]
    
    plot_locations_map(lons,lats,skaAntPosCsvFile=skaAntPosCsvFile)
    


