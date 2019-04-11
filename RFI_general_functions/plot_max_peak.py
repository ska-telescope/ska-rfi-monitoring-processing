# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 12:18:58 2019

@author: f.divruno
"""
import numpy as np
import matplotlib.pyplot as plt

#%%
# A class that will downsample the data and recompute when zoomed in a figure.
# The function plot_max_peak() should be used to plot very big amounts of data.
    
class DataDisplayDownsampler(object):
    def __init__(self, xdata, ydata):
        self.origYData = ydata
        self.origXData = xdata
        self.max_points = 1000
        self.delta = xdata[-1] - xdata[0]

    def downsample(self, xstart, xend):
        # get the points in the view range
        mask = (self.origXData >= xstart) & (self.origXData <= xend)
        # dilate the mask by one to catch the points just outside
        # of the view range to not truncate the line
        mask = np.convolve([1, 1], mask, mode='same').astype(bool)
        # sort out how many points to drop
        ratio = max(np.sum(mask) // self.max_points, 1)

        # mask data
        xdata = self.origXData[mask]
        ydata = self.origYData[mask]

        # downsample xdata
        xdata = xdata[::ratio]
        # calculate max peak for y data for every "ratio" number of samples.
        ydata = np.reshape(ydata,[len(ydata)//ratio,ratio])
        ydata = np.max(ydata,1)

        print("using {} of {} visible points".format(
            len(ydata), np.sum(mask)))

        return xdata, ydata

    def update(self, ax):
        # Update the line
        lims = ax.viewLim
        if np.abs(lims.width - self.delta) > 1e-8:
            self.delta = lims.width
            xstart, xend = lims.intervalx
            self.line.set_data(*self.downsample(xstart, xend))
            ax.figure.canvas.draw_idle()

def plot_max_peak(xdata,ydata):
    d = DataDisplayDownsampler(xdata, ydata)
    fig, ax = plt.subplots()

    # Hook up the line
#    d.line, = ax.plot(xdata, ydata)
    d.line = ax.plot(xdata, ydata) 
    ax.set_autoscale_on(False)  # Otherwise, infinite loop
    
    # Connect for changing the view limits
    ax.callbacks.connect('xlim_changed', d.update)
    ax.set_xlim(xdata.min(), xdata.max())
    plt.show()
    return d