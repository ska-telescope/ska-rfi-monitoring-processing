# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 10:43:40 2019

@author: f.divruno
"""

import visa
import matplotlib.pyplot as plt
import time
import numpy as np
from matplotlib.lines import Line2D
import matplotlib.animation as animation

class SpectAn(object):
    def __init__(self, ax,freq):
        self.ax = ax
        self.fdata = freq
        self.ydata = np.zeros(np.size(self.fdata))
        self.line = Line2D(self.fdata, self.ydata)
        self.ax.add_line(self.line)
        self.ax.set_ylim(-.1, 1.1)
        self.ax.set_xlim(self.fdata[0], self.fdata[-1])

    def update(self, y):
        self.ydata = y
        self.ax.set_ylim(np.min(self.ydata)-10, np.max(self.ydata)+20)
        #self.ax.figure.canvas.draw()
        self.line.set_data(self.fdata, self.ydata)
        return self.line,


def read_value(inst):
    while True:
        #inst.write("CLRW TRA;")
        #inst.query("DONE?;")
        #inst.write("MXMH TRA;")
        #inst.query("DONE?;")   
        #time.sleep(2)
        y = inst.query_ascii_values("TRA?;")
        inst.query("DONE?;")
        yield y


rm = visa.ResourceManager()
inst = rm.open_resource('GPIB0::18::INSTR')
print(inst.query("*ID?"))
inst.timeout = 100000
#inst.write("IP;")
#inst.query("DONE?;")

inst.write("TDF P;")

CF = 2450 #5785  #MHz
SPAN = 200 #MHz
REF = -40
ATT = 0
BW = 100 #kHz
Max_hold = 1
Single_sw = 0
Num_sw = 20

f = np.linspace(CF-(SPAN/2),CF+(SPAN/2),401)

inst.write("CF "+ str(CF)+ "Mz;")
inst.query("DONE?;")
inst.write("SP "+ str(SPAN)+ "MZ;")
inst.query("DONE?;")
inst.write("ATT "+ str(ATT)+ ";")
inst.query("DONE?;")
inst.write("RL "+ str(REF)+ ";")
inst.query("DONE?;")
inst.write("BW "+ str(BW)+ "Kz;")
inst.query("DONE?;")

if Max_hold:
     inst.write("MXMH TRA;")
     inst.query("DONE?;")

if Single_sw:
     inst.write("SNGLS;")
     inst.query("DONE?;")




fig, ax = plt.subplots()
SpectAn = SpectAn(ax,f)

# pass a generator in "emitter" to produce data for the update func
ani = animation.FuncAnimation(fig, SpectAn.update, read_value(inst), interval=5000,
                              blit=True)

plt.show()

