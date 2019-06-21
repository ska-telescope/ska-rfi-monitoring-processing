# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 00:04:14 2019
    Collection of functions 

@author: f.divruno
"""

import astropy.coordinates as Coord
import numpy as np

def Coord_to_nparray(coord):
    
    A = np.array([coord[0].si.value,coord[1].si.value,coord[2].si.value])
    
    return A