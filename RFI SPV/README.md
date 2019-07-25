# What is this repo?
- A framework for generation RFI signals and verification of their impact on complex systems such as the SKA  


# usage
```

Steps:
1) Copy an existing test case such as dmeSsrAircraftTestCase
2) Modify the following parameters to suit your situation (see below)
3) Add RFI emitters
4) Add telescopes
5) Run program

Parameters
----------
#file parameters
testCaseName = 'test1'
skaMidAntPosFileSpec = './skaMidAntPositions.csv'
randomSeed = 22.

#antenna pair to test
tstAnt1Key = 'SKA001'
tstAnt2Key = 'SKA005'

#antenna pointing az and el
antAzEl = dict(Elev=90*u.deg,Azimuth=0*u.deg)

#Receiver and temporal parameters
Band = 'B2'
Duration = 2.*ms
SamplingRate = 4*GHz # THis is the analog sampling rate

#ADC scaling
scaling = 'Correlator_opimized'

#Test configuration parameters
promptFlg = False #interactive mode prompts user at various processing steps
runFlg = True #can be used to skip over processing
saveFlg = True #results are saved if true
loadFlg = False #results are loaded if true
plot_signal = False #plot time series signal
plot_spectrum = False #plot spectrum
plot_corr = False   #plot correlation


TODO
add verification code and system performance metrics.
add main and parse parameters so program can be run from the command line
put parameters, emitters, antenna/receivers in a text file so the basic program structure can be maintained in one file.


```
