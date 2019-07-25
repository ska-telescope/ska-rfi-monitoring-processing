# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 21:43:25 2019
probability of detection of RFI

@author: f.divruno
"""


import numpy as np
import matplotlib.pyplot as plt

#%% 
puntos = 2000
#T1 = 160e-6 # RFI capture window
#T2 = 10 # RFI capture period

T1 = 100e-3 # RFI capture window
T2 = 1 # RFI capture period


#offset of RFI signal respect to the measurement period.
sig = T2
mu = 0
RFI_d = abs( sig*np.random.rand(puntos)+mu) #0 to T2 seconds, rectangular distribution.

#period of the RFI:
sig = 10 #seconds
mu = 2 #seconds
RFI_T = abs( sig*np.random.randn(puntos)+mu)

#T on of the RFI:
sig = 3 #seconds
mu = 0 #seconds
RFI_Ton = abs( sig*np.random.randn(puntos)+mu)


#number of pulses of the RFI:
sig = 5 #
mu = 3 #
RFI_M = np.round(abs( sig*np.random.randn(puntos))+mu)


Jmin = 0
Jmax = np.ceil((RFI_M*RFI_T+RFI_d+T2)/T2).astype(int)

R = np.zeros(puntos)

for i in range(puntos):
    print(i)
    J = np.linspace(0,Jmax[i],Jmax[i]+1)
    K = np.linspace(0,RFI_M[i].astype(int),RFI_M[i].astype(int)+1)
    [Jg,Kg] = np.meshgrid(J,K)
    result1 = (RFI_d[i]+Kg*RFI_T[i]) < (Jg*T2 + T1)
    result2 = (RFI_d[i]+Kg*RFI_T[i]+RFI_Ton[i]) > (J*T2)
    result = result1 * result2
    if np.sum(result[:].astype(int))>0:
        R[i] = 1

print('Total de detecciones = ' + str(sum(R)/puntos*100) +' %')   


#%% 

plt.subplot(2,2,1)
plt.hist(RFI_d,100)
plt.title('Histogram RFI delta')

plt.subplot(2,2,2)
plt.hist(RFI_Ton,100)
plt.title('Histogram RFI Ton')

plt.subplot(2,2,3)
plt.hist(RFI_T,100)
plt.title('Histogram RFI period')

plt.subplot(2,2,4)
plt.hist(RFI_M,100)
plt.title('Histogram RFI number of pulses')



