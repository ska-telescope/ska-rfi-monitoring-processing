"""
Created on Mon Jun 24 10:52:25 2019
    Reads a wav file with SDR IQ capture of FM stations located in :
    https://mega.nz/#F!3UUUnSiD!WLhWZ3ff4f4Pi7Ko_zcodQ
    
    Also generates IQ stream sampled at 2.4Msps to simulate a similar spectrum 
    sinusoids, this might be useful in an early stage to use a known signal.

@author: f.divruno
"""

#!/usr/bin/env python3


import wave
import numpy as np
import matplotlib.pyplot as plt


# ------------  PARAMETERS

N = 5000 #number of samples to read
nAverages = 10 # number of averages
folder = "C:\\Users\\F.Divruno\\Downloads\\" # change this to your folder.
filename = "17-22-08_89100kHz.wav"
CenterFrequency = 89100e3 # Centre freq of the recording is the number at the end of the filename.

# ------------




#Read an IQ recording of FM stations: 
wav_in = wave.open(folder+ filename, "r")
sampleFreq = 2.4e6 # sample freq of the SDR to acquire this signals

timeMax = N/sampleFreq # duration of the loaded signals
t = np.linspace(0,timeMax,N)

# Read the file
I = np.zeros(N)
Q = np.zeros(N)
for n in range(N):    
    aux = wav_in.readframes(1)  
    I[n] = aux[0]
    Q[n] = aux[1]


# Plot the spectrum of the recording
I_fft = np.fft.fftshift(np.fft.fft(I))
Q_fft = np.fft.fftshift(np.fft.fft(Q))
V = abs(I_fft-1j*Q_fft)

freq = np.fft.fftshift(np.fft.fftfreq(N,d=1/sampleFreq) + CenterFrequency)

plt.figure()
plt.subplot(2,1,1)
plt.plot(freq/1e6,20*np.log10(V))
plt.xlabel('MHz')
plt.ylabel('dB')
plt.title('Recording')



#test signal generated with tone signals
I = np.zeros(N)
Q = np.zeros(N)
foStation = np.array([88, 88.4, 88.6, 88.8 ,89.4, 89.6, 89.8, 90, 90.2])*1e6
numStations = len(foStation)
for k in range(numStations):
    fcent = (foStation[k] - CenterFrequency)/2
    for i in range(20):
        fc = fcent-10*1e3 + 1e3*i
        phase = np.random.random(1)*2*np.pi
        I += np.sin(2*np.pi*fc*t+phase)*np.sin(2*np.pi*fc*t)
        Q += np.sin(2*np.pi*fc*t+phase)*np.cos(2*np.pi*fc*t)


I_fft = np.fft.fftshift(np.fft.fft(I))
Q_fft = np.fft.fftshift(np.fft.fft(Q))
V = abs(I_fft-1j*Q_fft)

plt.subplot(2,1,2)
plt.plot(freq/1e6,20*np.log10(V),'g')
plt.xlabel('MHz')
plt.ylabel('dB')
plt.title('syntethized')


#%%Average the IQ recording of FM stations: 
wav_in.rewind()
V = np.zeros(N)
for k in range(nAverages):
    for n in range(N):    
        aux = wav_in.readframes(1)  
        I[n] = aux[0]
        Q[n] = aux[1]
    I_fft = np.fft.fftshift(np.fft.fft(I))
    Q_fft = np.fft.fftshift(np.fft.fft(Q))
    V += abs(I_fft-1j*Q_fft)

V /= nAverages


plt.figure()
plt.plot(freq/1e6,20*np.log10(V))
plt.title('Averaged %d times'%nAverages)
plt.xlabel('MHz')
plt.ylabel('dB')

