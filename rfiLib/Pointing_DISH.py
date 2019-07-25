# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 10:57:54 2019

    Ponting functions, starts from the elevation and azimuth coordinates and transforms them
    into the ECEF coordinates to use them for direction of arrival calculation.
    
@author: f.divruno
"""

#Rotation matrices
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

''' see https://gssc.esa.int/navipedia//index.php/Transformations_between_ECEF_and_ENU_coordinates
    for details on the transformation from the different reference frames.
'''

def ENU_to_ECEF(lat,lon):
    ''' returns the Transformation matrix to pass from the ENU
    reference frame to the ECEF frame.
    '''
    try:
        test1 = (lat.unit=='deg') & (lon.unit=='deg')
        if test1:
            sL = np.sin(lon)
            cL = np.cos(lon)
            sG = np.sin(lat)
            cG = np.cos(lat)
            
            T = np.array([[-sL,  -cL*sG, cL*cG],
                          [cL ,  -sL*sG, sL*cG],
                          [0  ,      cG,   sG]])

        else:
            raise Exception

    except:
        print('\nERROR ENU_to_ECEF:\nNeed to use deg units\n')
        raise SystemExit
    
    return T

def Pointing_to_ENU(elev,Az):
    '''transforms the elevation and azimuth parameters of pointing to
    the ENU coordinates in the ENU reference frame
    '''
    try:
        test1 = (elev.unit=='deg')
        test2 = (Az.unit=='deg')
        if test1 & test2:
            Penu = np.array([np.cos(elev)*np.cos(Az),
                             np.cos(elev)*np.sin(Az),
                             np.sin(elev)])
        else:
            raise Exception
    except:
        print('\nERROR Pointing_to_ENU:\nNeed to use deg units\n')
        raise SystemExit
       
    return Penu    


def ENU_to_latlon(P):
    '''
    transforms the ENU coordinates to elevation and azimuth 
    '''
    try:
        test1 = (np.size(P)==3)#(P.unit=='m')
        if test1:
            P_norm = np.linalg.norm(P)
            lat = np.arcsin(P[2]/P_norm)*180/np.pi
            lon = np.arctan2(P[1],P[0])*180/np.pi
        else:
            raise Exception
    except:
        print('ERROR ENU_to_latlon: Need to use m units')
        raise SystemExit
    return lat*u.deg,lon*u.deg


def Pointing_to_ECEF(elev,Az,lat,lon):
    '''transforms the elevation and azimuth parameters of pointing to
    the XYZ coordinates in the ECEF reference frame
    for repeated calculations is adviced to calculate the Transformation
    matrix from ENU to ECEF and then use it in a loop as in the built in test.
    Calling this function repeatedly will calculate the T matrix each time.
    '''
    try:
        test1 = (elev.unit=='deg') & (Az.unit=='deg') & (lat.unit=='deg') & (lon.unit=='deg')
        if test1:
            Penu = Pointing_to_ENU(elev,Az) # get the ENU coordinates           
            T = ENU_to_ECEF(lat,lon)    #Transformation matrix from ENU to ECEF
            Pecef = np.matmul(T,Penu)  #get the ECEF coordinates

        else:
            raise Exception
    except:
        print('\nERROR Pointing_to_ECEF:\nNeed to use deg units\n')
        raise SystemExit
            
    return Pecef



#built in test:
if __name__ == '__main__':

#    SKA_core = [5109.229, 2006.715,-3239.068]   #ECEF coordinates of the SKA MID site
#    S_norm = np.linalg.norm(SKA_core)           #magnitude
#    S = SKA_core/S_norm                         #ECEF versor
    
    #from Geodetic coordinates
    lat = -30.71278*u.deg      #Lat of the SKA-MID site
    lon = 21.44305*u.deg        #Long of the SKA-MID site
    h = 1000
    
    #ECEF :
#    Sx = np.cos(lat)*np.cos(lon)
#    Sy = np.cos(lat)*np.sin(lon)
#    Sz = np.sin(lat)
#    
    
    T = ENU_to_ECEF(lat,lon)    #Transformation matrix from ENU to ECEF
    
    Po = Pointing_to_ENU(lat,lon)   #ECEF versor of the SKA-MID site.
    

#    Plot of different pointing directions, the versor is drawn on the tip of the versor 
#    of the SKA-MID site to have a better visualization.
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot([0,Po[0]],[0,Po[1]],[0,Po[2]],'-b')
    
    for i in range(10):
        Az = i*36*u.deg
        for j in range (10):
            elev = 9*j*u.deg
            
            Penu = Pointing_to_ENU(elev,Az)
            P1 = np.matmul(T,Penu)
            
        #    ax.plot([S[0],S[0]+P1[0]],[S[1],S[1]+P1[1]],[S[2],S[2]+P1[2]],'-r')
            ax.plot([Po[0],Po[0]+P1[0]],[Po[1],Po[1]+P1[1]],[Po[2],Po[2]+P1[2]],label = 'Az= %.2f El=%.2f'%(Az.value,elev.value))
            
            cos_alpha = np.dot(P1,Po)
            print(np.arccos(cos_alpha)*180/np.pi)
        
#    ax.legend()
    
    