"""
The purpose of this python script is to plot scatter plot data as a contour plot
"""

import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from matplotlib import ticker


def plotContourFromScatter(x,y,ax,nbins=10,figsize=(6,8), norm=None,PPoutpath='./'):
    #norm could be LogNorm()
    norm=LogNorm()
    #DELETE x = np.random.normal(5,10,100000)
    #DELETE y = np.random.normal(5,10,100000)
    #fig = plt.figure(figsize=(6,8))
    #plt.rc('axes',linewidth=2)
    #plt.rc('lines',linewidth=2)
    #plt.rcParams['axes.linewidth']=2
    #plt.rc('font',weight='bold')
    delFig = plt.figure(num=65484)
    counts,ybins,xbins,image = plt.hist2d(x,y,bins=nbins,norm=norm,density=True)
    xcents = (xbins[:-1]+xbins[1:])/2.
    ycents = (ybins[:-1]+ybins[1:])/2.
    plt.close(65484)
    #norm=LogNorm()
    ax.contourf(xcents,ycents,counts,extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()], cmap='BuPu', norm=norm)#,locator=ticker.LogLocator())#'inferno')
    print('min: ' + str(np.min(counts)))
    print('max: ' + str(np.max(counts)))
    print('len xcents: ' + str(len(xcents)))
    print(xcents)
    #cmap='plasma' looked good
    #plt.show(block=False) 
    return ax

def plotKDEfromScatter(x,ax):
    #fig = plt.figure(figsize=(6,6))
    kernel = 'epanechnikov'
    kde = KernelDensity(kernel=kernel, bandwidth=(np.max(x)-np.min(x))/50.).fit(np.asarray(x).reshape(-1, 1))
    X_plot = np.linspace(np.min(x), np.max(x), 1000)[:, np.newaxis] #the range of values to plot over
    log_dens = kde.score_samples(X_plot) # the plotted values
    ax.plot(X_plot[:, 0], np.exp(log_dens), '-', color='black')
            #label="kernel = '{0}'".format(kernel))
    #plt.show(block=False)
    return ax

