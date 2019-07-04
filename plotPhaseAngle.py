"""
The purpose of this script is to take a spacecraft TLE
and correlate a specific ground observation's photometric and astrometric
uncertainties to quantify the spacecraft size
"""

import numpy as np
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
import matplotlib.pyplot as plt
from astropy import constants as const
import copy

# For a given TLE, try running 
# sgp4 using from sgp4.io import twoline2rv
#and pyorbital using from pyorbital import tlefile

def extractTLEdata(TLEfilepath,TLEfilename):
    """ Extract data from TLEs
    Args:
        TLEfilepath (string) - The filepath to where TLEfilename is stored
        TLEfilename (string) - The TLE file name
    """
    #### Open space-track TLE file ###################################
    import os
    #DELETE with open("/home/dean/Documents/AFRL2019/tle-new.txt", 'r') as f:
    #DELETE with open("/home/dean/Documents/AFRL2019/3le.txt", 'r') as f:
    with open(os.path.join(TLEfilepath,TLEfilename), 'r') as f:
        # read all TLE from file and remove final line caused by split
        lines = f.read().split('\n')[:-1] 

    #### Parsing space-track TLE data ################################
    assert np.mod(len(lines),3) == 0, 'the number of lines for each TLE is not 3'
    numTLE = len(lines)/3 # there should be 1 TLE for every 3 lines

    satData = {} # stores all parsed satellite data
    satNums = list() # stores all satellite identification numbers
    for i in np.arange(numTLE, dtype='int'):
        line = lines[i*3]
        assert line[0] == '0', 'Error: Not a Name Line'
        #The first character on every name line must be 0    
        satName = line.replace(" ","")
        TLE1 = lines[i*3+1] # first line of TLE
        TLE2 = lines[i*3+2] # second line of TLE

        satOBJ = twoline2rv(TLE1, TLE2, wgs72)
        satNums.append(satOBJ.satnum)
        satData[satOBJ.satnum] = {
        "satName":satName,
        "satelliteOBJ":satOBJ,
        "TLE0":line,
        "TLE1":TLE1,
        "TLE2":TLE2
        }
    return satData, lines, numTLE

#### Open space-track TLE file #########################################
TLEfilepath = "/home/dean/Documents/AFRL2019"
TLEfilename = "3le.txt"
satData, lines, numTLE = extractTLEdata(TLEfilepath,TLEfilename)
print('Done Loading TLE Data')
del TLEfilepath, TLEfilename
########################################################################

##### Extracting individual TLE orbital parameters #####################
no = [satData[key]['satelliteOBJ'].no for key in satData.keys()]#satOBJ.no - mean motion at epoch in rad/min
ecco = [satData[key]['satelliteOBJ'].ecco for key in satData.keys()]#satOBJ.ecco - eccentricity
inclo = [satData[key]['satelliteOBJ'].inclo for key in satData.keys()]#satOBJ.inclo - mean inclination in rad
mo = [satData[key]['satelliteOBJ'].mo for key in satData.keys()]#satOBJ.mo - mean Mean Anomaly in rad
argpo = [satData[key]['satelliteOBJ'].argpo for key in satData.keys()]#satOBJ.argpo - mean argument of perigee in rad
nodeo = [satData[key]['satelliteOBJ'].nodeo for key in satData.keys()]#satOBJ.nodeo - mean longitude of the ascending node 
#       (right ascension of the ascending node) in rad
a = [satData[key]['satelliteOBJ'].a for key in satData.keys()]#satOBJ.a - semimajor-axis in Earth radii
alta = [satData[key]['satelliteOBJ'].alta for key in satData.keys()]#satOBJ.alta - altitude at apoapsis in Earth radii
altp = [satData[key]['satelliteOBJ'].altp for key in satData.keys()]#satOBJ.altp - altitude at periapsis in Earth radii
#epocyr = [satData[key]['satelliteOBJ'].no for key in satData.keys()]#satOBJ.epocyr - year of the epoch XXXX
jdsatepoch = [satData[key]['satelliteOBJ'].jdsatepoch for key in satData.keys()]#satOBJ.jdsatepoch - time since J2000 in julian days
#epoch = [satData[key]['satelliteOBJ'].no for key in satData.keys()]#satOBJ.epoch - time the snapshot was taken
print('Done Extracting Individual TLE Parameters')
del satData, lines, numTLE
########################################################################

#### Calculate Eccentric Anomaly, E ####################################
from EXOSIMS.util.eccanom import eccanom
E = [eccanom(mo[i], ecco[i])[0] for i in np.arange(len(mo))] #eccentric anomaly
########################################################################

#### Plot Scatter Matrix ###############################################
from pandas.plotting import scatter_matrix
import pandas as pd 

#Format data to turn into pandas frame
pdData = {'Mean\nMotion\n(rad/min)':no, 'Eccentricity':ecco, 'Inclination\n(rad)':inclo,\
    'Mean\nAnomaly\n(rad)':mo, 'Arg. of\nPerigee\n(rad)':argpo, 'Longitude\nAscending\nNode\n(rad)':nodeo,\
    'Semi-major\nAxis\n(Earth Radii)':a,'Apoapsis\nAltitude\n(Earth Radii)':alta, 'Periapsis\nAtlitude\n(Earth radii)':altp,\
    'Time Since\nEpoch (JD)':jdsatepoch, 'Eccentric\nAnomaly\n(rad)':E}
del no, ecco, inclo, mo, argpo, nodeo, a, alta, altp, jdsatepoch, E

#### Plot Orbit Dist Scatter Matrix #########################################
# def plotOrbitDist_scatterMatrix(pdData):
# """ Plots a scatter matrix of Keplerian Orbital Parameter Data
# Args:
#     pdData (pandas dataframe object)
# """
df = pd.DataFrame(data=pdData)
fignum = 651686214
plt.close(fignum)
fig = plt.figure(num=fignum,figsize=(12,12))
plt.rc('axes',linewidth=2)
plt.rc('lines',linewidth=2)
plt.rcParams['axes.linewidth']=2
plt.rc('font',weight='bold')
ax = plt.axes()
ax2 = scatter_matrix(df, alpha=0.01, diagonal='kde', ax=ax, **{'color':'black'})#, **kwds)
for ax_sub1 in ax2:
    for ax_sub2 in ax_sub1:
        label = ax_sub2.get_ylabel()
        ax_sub2.set_ylabel(label,rotation=0, labelpad=40, weight='bold')
        if 'Motion' in label: # the mean motion label has bad tickmarkers 0.050000000001 or something like that
            tmplabels = ax_sub2.get_ymajorticklabels()
            if not tmplabels[0].get_text() == '': #cant be empty for float() to work
                for i in np.arange(len(tmplabels)):
                    txt = tmplabels[i].get_text()
                    tmplabels[i].set_text("{:.3f}".format(np.round(float(txt),decimals=3)))
                    del txt
                ax_sub2.set_yticklabels(tmplabels)
            del tmplabels
        label2 = ax_sub2.get_xlabel()
        ax_sub2.set_xlabel(label2, weight='bold')
        del label2
plt.show(block=False)
print('Done Plotting Scatter Matrix')
plt.pause(2.)
del df

#OK### Plot Eccen. vs SMA relation based on R_earth Min
#EXCLUSION REGION IS TOP LEFT
SMA = np.linspace(start=np.min(pdData['Semi-major\nAxis\n(Earth Radii)']),
        stop=np.max(pdData['Semi-major\nAxis\n(Earth Radii)']), num=30) # Semi-Major Axis
R_earth = 1. #In earth radii
minAlt = 0.#LEO min alt is technically earth surface 400./6356. #assumes 400km is min LEO alt, Earth Radius is 6356 km from Wikipedia
R = R_earth + minAlt
ECCEN = 1. - R/SMA
ax2[1][6].plot(SMA,ECCEN,color='b',alpha=0.4)
plt.show(block=False)
#OK### Plot SMA vs Eccen. relation based on R_earth Min
#EXCLUSION REGION IS BOTTOM RIGHT
ECCEN = np.linspace(start=np.min(pdData['Eccentricity']),
        stop=np.max(pdData['Eccentricity']), num=30) # Eccentricity
R_earth = 1. #In earth radii
minAlt = 0.#LEO min alt is technically earth surface 400./6356. #assumes 400km is min LEO alt, Earth Radius is 6356 km from Wikipedia
R = R_earth + minAlt
SMA = R/(1.-ECCEN)
ax2[6][1].plot(ECCEN,SMA,color='b',alpha=0.4)
plt.show(block=False)
#OK### Plot Eccen. vs Apoapsis relation based on R_earth Min
#EXCLUSION REGION IS TOP LEFT
APOAPSIS = np.linspace(start=np.min(pdData['Apoapsis\nAltitude\n(Earth Radii)']),
        stop=np.max(pdData['Apoapsis\nAltitude\n(Earth Radii)']), num=30) # Apoapsis
R_earth = 1. #In earth radii
minAlt = 0.#LEO min alt is technically earth surface 400./6356. #assumes 400km is min LEO alt, Earth Radius is 6356 km from Wikipedia
R = R_earth + minAlt
ECCEN = 2.*APOAPSIS/(R+APOAPSIS) - 1.
ax2[1][7].plot(APOAPSIS-R_earth,ECCEN,color='b',alpha=0.4)
plt.show(block=False)
#OK### Plot Apoapsis vs Eccen. relation based on R_earth Min
#EXCLUSION REGION IS BOTTOM RIGHT
ECCEN = np.linspace(start=np.min(pdData['Eccentricity']),
        stop=np.max(pdData['Eccentricity']), num=30) # Eccentricity
R_earth = 1. #In earth radii
minAlt = 0.#LEO min alt is technically earth surface 400./6356. #assumes 400km is min LEO alt, Earth Radius is 6356 km from Wikipedia
R = R_earth + minAlt
APOAPSIS = R*(ECCEN+1.)/(1.-ECCEN)
ax2[7][1].plot(ECCEN,APOAPSIS-R_earth,color='b',alpha=0.4)
plt.show(block=False)
#OK### Plot APOAPSIS vs SMA
SMA = np.linspace(start=np.min(pdData['Semi-major\nAxis\n(Earth Radii)']),
        stop=np.max(pdData['Semi-major\nAxis\n(Earth Radii)']), num=30) # Semi-Major Axis
APOAPSIS = SMA-R
ax2[7][6].plot(SMA,APOAPSIS,color='r',alpha=0.4)
plt.show(block=False)
#OK### Plot SMA vs APOAPSIS
APOAPSIS = np.linspace(start=np.min(pdData['Apoapsis\nAltitude\n(Earth Radii)']),
        stop=np.max(pdData['Apoapsis\nAltitude\n(Earth Radii)']), num=30) # Apoapsis
SMA = APOAPSIS+R
ax2[6][7].plot(APOAPSIS,SMA,color='r',alpha=0.4)
plt.show(block=False)
#OK### Plot PERIAPSIS min on ECCEN vs PERIAPSIS plot, since periapsis can't be smaller than the earth radius
R_earth = 1. #In earth radii
minAlt = 0.#LEO min alt is technically earth surface 400./6356. #assumes 400km is min LEO alt, Earth Radius is 6356 km from Wikipedia
R = R_earth + minAlt
ax2[1][8].plot([R-R,R-R],[ax2[1][8].get_xlim()[0],ax2[1][8].get_xlim()[1]],color='r',alpha=0.4)
plt.show(block=False)
#OK### Plot PERIAPSIS min on PERIAPSIS vs ECCEN, since periapsis can't be smaller than the earth radius
R_earth = 1. #In earth radii
minAlt = 0.#LEO min alt is technically earth surface 400./6356. #assumes 400km is min LEO alt, Earth Radius is 6356 km from Wikipedia
R = R_earth + minAlt
ECCEN = np.linspace(start=np.min(pdData['Eccentricity']),
        stop=np.max(pdData['Eccentricity']), num=30) # Eccentricity
ax2[8][1].plot(ECCEN,R-R + np.zeros(len(ECCEN)),color='r',alpha=0.4)
plt.show(block=False)
#OK### Plot Minimum APOAPSIS vs PERIAPSIS, since APOAPSIS must be larger than PERIAPSIS #BOTTOM RIGHT EXCLUSION ZONE
PERIAPSIS = np.linspace(start=np.min(pdData['Periapsis\nAtlitude\n(Earth radii)']),
        stop=np.max(pdData['Periapsis\nAtlitude\n(Earth radii)']), num=30) # Periapsis
APOAPSIS = PERIAPSIS
ax2[7][8].plot(PERIAPSIS,APOAPSIS,color='r',alpha=0.4)
plt.show(block=False)
#OK### Plot MAXIMUM PERIAPSIS vs APOAPSIS, since APOAPSIS must be larger than PERIAPSIS #TOP LEFT EXCLUSION ZONE
APOAPSIS = np.linspace(start=np.min(pdData['Apoapsis\nAltitude\n(Earth Radii)']),
        stop=np.max(pdData['Apoapsis\nAltitude\n(Earth Radii)']), num=30) # Apoapsis
PERIAPSIS = APOAPSIS
ax2[8][7].plot(APOAPSIS,PERIAPSIS,color='r',alpha=0.4)
plt.show(block=False)
#OK### Plot Mean Motion vs SMA
SMA = np.linspace(start=np.min(pdData['Semi-major\nAxis\n(Earth Radii)']),
        stop=np.max(pdData['Semi-major\nAxis\n(Earth Radii)']), num=30) # Semi-Major Axis
mu = const.GM_earth
MEANMOTION = np.sqrt(mu.value*(60.*60.)/(const.R_earth.value**3.)/SMA**3.)
ax2[0][6].plot(SMA,MEANMOTION,color='r',alpha=0.4)
plt.show(block=False)
#OK### Plot SMA vs Mean Motion
MEANMOTION = np.linspace(start=np.min(pdData['Mean\nMotion\n(rad/min)']),
        stop=np.max(pdData['Mean\nMotion\n(rad/min)']), num=30) # Mean Motion
mu = const.GM_earth
MEANMOTION = np.sqrt(mu.value*(60.*60.)/(const.R_earth.value**3.)/SMA**3.)
SMA = (mu.value*(60.*60.)/(const.R_earth.value**3.)/MEANMOTION**2.)**(1./3.)
ax2[6][0].plot(MEANMOTION,SMA,color='r',alpha=0.4)
plt.show(block=False)
#OK### Plot Mean Motion vs APOAPSIS
R_earth = 1.
APOAPSIS = np.linspace(start=np.min(pdData['Apoapsis\nAltitude\n(Earth Radii)']),
        stop=np.max(pdData['Apoapsis\nAltitude\n(Earth Radii)']), num=30) # Apoapsis
MEANMOTION = np.sqrt(mu.value*(60.*60.)/(const.R_earth.value**3.)/(APOAPSIS)**3.)
ax2[0][7].plot(APOAPSIS-R_earth,MEANMOTION,color='r',alpha=0.4) #should be a lower bound
MEANMOTION = np.sqrt(8.*mu.value*(60.*60.)/(const.R_earth.value**3.)/(APOAPSIS)**3.) #should be an upper limit
ax2[0][7].plot(APOAPSIS-R_earth,MEANMOTION,color='r',alpha=0.4)
plt.show(block=False)
#NOT TOTALLY SURE ABOUT THE X-AXIS VALUES. NEED TO DOUBLE CHECK WHETHER APOAPSIS HERE IS CORRECT
#OK### Plot Mean Motion vs PERIAPSIS
R_earth = 1.
PERIAPSIS = np.linspace(start=np.min(pdData['Periapsis\nAtlitude\n(Earth radii)']),
        stop=np.max(pdData['Periapsis\nAtlitude\n(Earth radii)']), num=30) # Periapsis
ax2[0][8].plot(PERIAPSIS-R_earth,0. + np.zeros(len(PERIAPSIS)),color='r',alpha=0.4) # since 1-e_max = 0 
MEANMOTION = np.sqrt(mu.value*(60.*60.)/(const.R_earth.value**3.)/(PERIAPSIS)**3.)
ax2[0][8].plot(PERIAPSIS-R_earth,MEANMOTION,color='r',alpha=0.4) # since 1-e_min = 0 
plt.show(block=False)

print('Done plotting Scatter Matrix Limits')
plt.pause(2.)
del SMA, MEANMOTION, APOAPSIS, PERIAPSIS, ECCEN, R, minAlt, R_earth, mu

#### Save Scatter Data
import datetime
import re
import os
# Save to a File
plotBOOL = False
if plotBOOL==True:
    PPoutpath = '/home/dean/Documents/AFRL2019'
    folder = PPoutpath
    date = str(datetime.datetime.now())
    date = ''.join(c + '_' for c in re.split('-|:| ',date)[0:-1])#Removes seconds from date
    fname = 'SpaceObjectOrbitalParameters_' + folder.split('/')[-1] + '_' + date
    plt.savefig(os.path.join(PPoutpath, fname + '.png'), format='png', dpi=200)
    plt.savefig(os.path.join(PPoutpath, fname + '.svg'))
    #plt.savefig(os.path.join(PPoutpath, fname + '.eps'), format='eps', dpi=500) #doesnt plot transparent stuff
    plt.savefig(os.path.join(PPoutpath, fname + '.pdf'), format='pdf', dpi=200)
    print('Done Saving Scatter Matrix Figure')
    del PPoutpath, folder, data, fname
else:
    print('Skipping Saving Scatter Matrix Figure')
###########################################################################

#comment-out to improve speed
#plotOrbitDist_scatterMatrix(pdData)
#############################################################################

#### Plot Spacecraft Parameters in "contour matrix" #######################
#### Generate Joint PDFs
from plotContourFromScatter import plotContourFromScatter
from plotContourFromScatter import plotKDEfromScatter
from scipy import interpolate

def get_pdDataLimits(pdData,key1,key2):
    """
    Consolidates pdData limit calculation and exceptions
    By default, data limit is the maximal or minimal limits of the provided data
    Otherwise, it is a known limit (like 0-2pi or 0-1)
    """
    xaxismin = np.min(pdData[key1])
    xaxismax = np.max(pdData[key1])
    yaxismin = np.min(pdData[key2])
    yaxismax = np.max(pdData[key2])
    #0 to 1 limits
    if key1 == 'Eccentricity':
        xaxismin = 0.
        xaxismax = 1.
    elif key2 == 'Eccentricity':
        yaxismin = 0.
        yaxismax = 1.
    #0 to 2pi limits
    if key1 in ['Mean\nAnomaly\n(rad)', 'Arg. of\nPerigee\n(rad)', 'Longitude\nAscending\nNode\n(rad)',\
        'Eccentric\nAnomaly\n(rad)']:
        xaxismin = 0.
        xaxismax = 2.*np.pi
    elif key2 in ['Mean\nAnomaly\n(rad)', 'Arg. of\nPerigee\n(rad)', 'Longitude\nAscending\nNode\n(rad)',\
        'Eccentric\nAnomaly\n(rad)']:
        yaxismin = 0.
        yaxismax = 2.*np.pi
    #0 to pi limits
    if key1 == 'Inclination\n(rad)':
        xaxismin = 0.
        xaxismax = np.pi
    elif key2 == 'Inclination\n(rad)':
        yaxismin = 0.
        yaxismax = np.pi
    return xaxismin, xaxismax, yaxismin, yaxismax

# Create Lists to Hold PDF and vectorized PDFS
mDFL = {}
for key1 in pdData.keys():
    mDFL[key1] = {}
    for key2 in pdData.keys():
        mDFL[key1][key2] = {}

#Create all axis options
fignum=516841
plt.close(fignum)
plt.rc('axes',linewidth=2)
plt.rc('lines',linewidth=2)
plt.rcParams['axes.linewidth']=2
plt.rc('font',weight='bold')
outAX = plt.subplots(len(pdData.keys()),len(pdData.keys()), num=fignum, figsize=(12,12))
outAX[0].subplots_adjust(bottom=0.15)
#outAX[0] is the figure object
for ikey1 in np.arange(len(pdData.keys())):#ROWS
    for ikey2 in np.arange(len(pdData.keys())):#COLUMNS
        key1 = list(pdData.keys())[ikey1]
        key2 = list(pdData.keys())[ikey2]
        if key1 == key2:
            outAX[1][ikey1][ikey2] = plotKDEfromScatter(pdData[key1], ax=outAX[1][ikey1][ikey1])
        else:
            outAX[1][ikey1][ikey2] = plotContourFromScatter(pdData[key1],pdData[key2],ax=outAX[1][ikey1][ikey2],nbins=11,figsize=(6,8))
            #### Create Joint Probability Distributions
            xaxismin, xaxismax, yaxismin, yaxismax = get_pdDataLimits(pdData,key2,key1)
            bins=134 # this is calculated by sqrt(17966)
            h, xedges, yedges = np.histogram2d(pdData[key2], pdData[key1], bins=bins, range=([xaxismin,xaxismax],[yaxismin,yaxismax]),density=False)
            h = np.pad(h,0,mode='constant')
            xcent = 0.5*(xedges[1:]+xedges[:-1])
            ycent = 0.5*(yedges[1:]+yedges[:-1])
            xnew = np.hstack((xaxismin,xcent,xaxismax))
            ynew = np.hstack((yaxismin,ycent,yaxismax))
            mDFL[key1][key2]['Hdist'] = h #stores NON-NORMALIZED 2d histogram
            mDFL[key1][key2]['xnew'] = xnew #stores xnew
            mDFL[key1][key2]['ynew'] = ynew #stores ynew
            mDFL[key1][key2]['xedges'] = xedges #stores xedges
            mDFL[key1][key2]['yedges'] = yedges #stores yedges
            mDFL[key1][key2]['xcent'] = xcent #stores xcent
            mDFL[key1][key2]['ycent'] = ycent #stores ycent
            hNORM, xedges, yedges = np.histogram2d(pdData[key2], pdData[key1], bins=bins, range=([xaxismin,xaxismax],[yaxismin,yaxismax]),density=True)
            #hNORM = np.pad(hNORM,0.,mode='constant')
            #normed is the probability density normalized by bin area. in HdistNorm, I multiply by bin areas
            mDFL[key1][key2]['HdistNORM'] = hNORM #stores 2d histogram
            #NOTE: I lightly stretch the axis in order to cover the entire 0 to 2pi or 0 to 1 space 
            normalizing_constant = interpolate.RectBivariateSpline(np.linspace(start=xedges[0],stop=xedges[-1],num=bins,endpoint=True),\
                np.linspace(start=yedges[0],stop=yedges[-1],num=bins,endpoint=True), hNORM, bbox=[xaxismin, xaxismax, yaxismin, yaxismax]).integral(-np.inf,np.inf,-np.inf,np.inf)
            EVPOCpdf = interpolate.RectBivariateSpline(np.linspace(start=xedges[0],stop=xedges[-1],num=bins,endpoint=True),\
                np.linspace(start=yedges[0],stop=yedges[-1],num=bins,endpoint=True), hNORM/normalizing_constant, bbox=[xaxismin, xaxismax, yaxismin, yaxismax])#*((xaxismax - xaxismin)*(yaxismax - yaxismin)))##/np.sum(h))
            #NOTE: Using kx,ky,s in RBS is a bad idea
            mDFL[key1][key2]['EVPOCpdf'] = copy.deepcopy(EVPOCpdf) #stores PDF
            mDFL[key1][key2]['vectorizedEVPOCpdf'] = np.vectorize(copy.deepcopy(EVPOCpdf))
            #DELETE if key2 == 'Eccentricity' and key1 == 'Eccentric\nAnomaly\n(rad)':
            #     print(saltyburrito)
            #keep EVPOC = np.vectorize(EVPOCpdf.integral, otypes=[np.float64])
            #keep mDFL[key1][key2]['EVPOCpdf'] = EVPOC #stores vectorized PDF
            #################################################################

        #Turns off all y-axis keys except first column
        if not key2 == list(pdData.keys())[0]: # the y-axis keys
            outAX[1][ikey1][ikey2].get_yaxis().set_visible(False)
        else:
            outAX[1][ikey1][ikey2].set_ylabel('\n' + list(pdData.keys())[ikey1],weight='bold',rotation=0,labelpad=40) #set ylabel
            #NEED ADJUSTMENT IF KEY2=0, USE Y AXIS LABELS OF PLOT TO RIGHT
            if key1 == list(pdData.keys())[0]: #this is the first key in the list
                tmplabels = outAX[1][ikey1][ikey2+1].get_ymajorticklabels()
                if not tmplabels[0].get_text() == '': #cant be empty for float() to work
                    for i in np.arange(len(tmplabels)):
                        txt = tmplabels[i].get_text()
                        tmplabels[i].set_text("{:.3f}".format(np.round(float(txt),decimals=3)))
                        del txt
                    outAX[1][ikey1][ikey2].set_yticklabels(tmplabels)
                del tmplabels
        #Turns off x-axis labels for all plots except bottom row
        if not key1 == list(pdData.keys())[-1]: # the x-axis keys
            outAX[1][ikey1][ikey2].get_xaxis().set_visible(False)
        else:
            outAX[1][ikey1][ikey2].tick_params(axis='x', rotation=90)
            outAX[1][ikey1][ikey2].set_xlabel(list(pdData.keys())[ikey2],weight='bold') #set xlabel
            
plt.subplots_adjust(wspace=0, hspace=0)
plt.show(block=False)
print('done creating all contour plots')
plt.pause(2.)
"""
mDFL contains fields key1 - x parameter, key2 - y parameter:
xedges - xedges of the bins from the 2d histogram
xnew - [xmin,xcent,xmax] of bins
xcent - bin centers of edges
yedges - xedges of the bins from the 2d histogram
ynew - [ymin,ycent,ymax] of bins
ycent - bin centers of edges
Hdist - the 2d histogram of quantity of spacecraft in each bin (sums to the number of SC ~17966)
HdistNORM - the 2d histogram of probability a spacecraft exists in a bin per bin area. (to get #sc, multiply by #sc and total area)
EVPOCpdf - the rectilinear bivariate spline (RBS) interpolant over the 2d histogram (can be queried at x,y)
vectorizedEVPOCpdf - the vectorized RBS interpolant which can be queries over arrays (x0...xn,y0...yn)
"""
############################################################################################################
#plt.close()

#OK### Calculate PDF and InverseTransform Sampler From data #################################################################################
from scipy.stats.kde import gaussian_kde
from numpy import linspace
from EXOSIMS.util.InverseTransformSampler import InverseTransformSampler

kdes = list()
sampler = list()
sampledVals = list()
for i in np.arange(len(pdData.keys())):
    key1 = list(pdData.keys())[i]
    kdes.append(gaussian_kde(pdData[key1]))
    #we now have a list of KDEs for each of the keplerian orbital parameters
    ########################################################################

    xmin = np.min(pdData[key1])#ax2[i][i].get_xlim()[0]
    xmax = np.max(pdData[key1])#ax2[i][i].get_xlim()[1]
    xedges = np.linspace(start=xmin,stop=xmax,num=100)
    cdf = np.vectorize(lambda x: kdes[i].integrate_box_1d(-np.inf, x))

    pdfvals = kdes[i].pdf(xedges)
    cdfvals = cdf(xedges)

    nsamples=50000
    sampler.append(InverseTransformSampler(kdes[i].pdf, xmin, xmax)) # contains the Inverse transform sampler for each parameter
    sampledVals.append(sampler[i](nsamples)) # contains the list of sampled values for each parameter
    #We now have a sample and samples for each parameter
print('Done Calculating Single-Variable Inverse Transform Sampler')
####################################################################################################################

#### Plot all Recilinear bivariate splines #########################################################################
#Create all axis options
fignum=516842
plt.close(fignum)
plt.rc('axes',linewidth=2)
plt.rc('lines',linewidth=2)
plt.rcParams['axes.linewidth']=2
plt.rc('font',weight='bold')
outAX2 = plt.subplots(len(pdData.keys()),len(pdData.keys()), num=fignum, figsize=(12,12))
outAX2[0].subplots_adjust(bottom=0.15)
#outAX2[0] is the figure object
for ikey1 in np.arange(len(pdData.keys())):
    for ikey2 in np.arange(len(pdData.keys())):
        key1 = list(pdData.keys())[ikey1]
        key2 = list(pdData.keys())[ikey2]
        if key1 == key2:
            outAX2[1][ikey1][ikey2] = plotKDEfromScatter(pdData[key1], ax=outAX2[1][ikey1][ikey1])
        else:
            # yaxismin = np.min(pdData[key1])
            # yaxismax = np.max(pdData[key1])
            # xaxismin = np.min(pdData[key2])
            # xaxismax = np.max(pdData[key2])
            # x2 = np.linspace(start=xaxismin,stop=xaxismax,num=100) #these are actually the Y-axis values
            # y2 = np.linspace(start=yaxismin,stop=yaxismax,num=100) # these are actually the X-axis values
            x2 = mDFL[key1][key2]['xcent']
            y2 = mDFL[key1][key2]['ycent']
            xgrid, ygrid = np.meshgrid(x2,y2) 
            Z2 = mDFL[key1][key2]['EVPOCpdf'](x2,y2)
            outAX2[1][ikey1][ikey2].contourf(xgrid,ygrid,Z2.T, levels=100, cmap='bwr')

        if not key2 == list(pdData.keys())[0]: # the y-axis keys
            outAX2[1][ikey1][ikey2].get_yaxis().set_visible(False)
        else:
            outAX2[1][ikey1][ikey2].set_ylabel('\n' + list(pdData.keys())[ikey1],weight='bold',rotation=0,labelpad=40) #set ylabel

            #NEED ADJUSTMENT IF KEY2=0, USE Y AXIS LABELS OF PLOT TO RIGHT
        if not key1 == list(pdData.keys())[-1]: # the x-axis keys
            outAX2[1][ikey1][ikey2].get_xaxis().set_visible(False)
        else:
            outAX2[1][ikey1][ikey2].tick_params(axis='x', rotation=90)
            outAX2[1][ikey1][ikey2].set_xlabel(list(pdData.keys())[ikey2],weight='bold') #set xlabel
            
plt.subplots_adjust(wspace=0, hspace=0)
plt.show(block=False)
print('done creating all RectilinearBivariateSpline plots')
plt.pause(2.)
#########################################################################################################################
#EXAMPLE FOR TESTING THE IMMEDIATELY ABOVE CODE DELETE WHEN DONE
# key1 = list(pdData.keys())[0]
# key2 = list(pdData.keys())[6]#[3]#[6]
# yaxismin = np.min(pdData[key1])
# yaxismax = np.max(pdData[key1])
# xaxismin = np.min(pdData[key2])
# xaxismax = np.max(pdData[key2])
# x2 = np.linspace(start=xaxismin,stop=xaxismax,num=100) #these are actually the Y-axis values
# y2 = np.linspace(start=yaxismin,stop=yaxismax,num=100) # these are actually the X-axis values
# xgrid, ygrid = np.meshgrid(x2,y2)
# Z2 = mDFL[key1][key2]['EVPOCpdf'](x2,y2)
# plt.close(24516843155)
# fig1213 = plt.figure(24516843155)
# plt.contourf(xgrid,ygrid,Z2.T, levels=100, cmap='bwr')#cmap='Greys')
# plt.show(block=False)
####

#DELETE THIS, THE ABOVE CODE ALREADY CALCULATES EVPOCPDF
# #### Calculate 2d distribution for Eccentric Anomaly E: Attempt 1 #############################################################
# Eind = [i for i in np.arange(len(list(pdData.keys()))) if list(pdData.keys())[i] == 'Eccentric\nAnomaly\n(rad)'][0] #assumes only 1 found
# eccind = [i for i in np.arange(len(list(pdData.keys()))) if list(pdData.keys())[i] == 'Eccentricity'][0] #assumes only 1 found

# Es = pdData[list(pdData.keys())[Eind]] # these are the calculated Eccentric Anomalies from the SC TLE data
# eccs = pdData[list(pdData.keys())[eccind]] # these are the calculated eccentricities from the SC TLE data

# xmin = np.min(Es)
# xmax = np.max(Es)
# ymin = np.min(eccs)
# ymax = np.max(eccs)
# Eshist = np.histogram(Es)
# eccshist = np.histogram(eccs)
# h, xedges, yedges = np.histogram2d(Es, eccs, bins=1000, range=([xmin,xmax],[ymin,ymax]),density=True)
# xcent = 0.5*(xedges[1:]+xedges[:-1])
# ycent = 0.5*(yedges[1:]+yedges[:-1])
# xnew = np.hstack((xmin,xcent,xmax))
# ynew = np.hstack((ymin,ycent,ymax))
# EeccPDF = h.T
# EeccPDF = np.pad(EeccPDF,1,mode='constant')
# EVPOCpdf = interpolate.RectBivariateSpline(xnew, ynew, EeccPDF)
# EVPOC = np.vectorize(EVPOCpdf.integral, otypes=[np.float64])
# levels = np.linspace(start=h.min(),stop=h.max(),num=100,endpoint=True)

# plt.close(1242511)
# fig = plt.figure(num=1242511)
# #plt.contourf(xnew[1:-1],ynew[1:-1],h)
# plt.scatter(Es,eccs,alpha=0.1)
# plt.show(block=False)

# plt.close(1242512)
# fig = plt.figure(num=1242512)
# plt.contourf(xnew[1:-1],ynew[1:-1],h.T, levels=levels)
# plt.show(block=False)


# #sampler[Eind]
# # Es = sampledVals[Eind] #the sampled values


# # xmin = ax2[1][6].get_xlim()[0]
# # xmax = ax2[1][6].get_xlim()[1]
# # xedges = np.linspace(start=xmin,stop=xmax,num=100)
# # xcent = 0.5*(xedges[1:]+xedges[:-1])
# # ymin = ax2[1][6].get_ylim()[0]
# # ymax = ax2[1][6].get_ylim()[1]
# # yedges = np.linspace(start=ymin,stop=ymax,num=100)
# # ycent = 0.5*(yedges[1:]+yedges[:-1])
# # xnew = np.hstack((0.0,xcent,xmax))
# # ynew = np.hstack((ymin,ycent,ymax))

# #NEED TO CALCULATE CPDF
# # Cpdf = np.pad(Cpdf,1,mode='constant')
# # #save interpolant to object
# # self.Cpdf = Cpdf
# # self.EVPOCpdf = interpolate.RectBivariateSpline(xnew, ynew, Cpdf.T)
# ####################################################################################################################

#### CALCULATE THE DISTRIBUTION OF MEAN ANOMALY




#### Calculate 2d distribution for Eccentric Anomaly E and ecentricity: Attempt 2 ##################################
#This generates a random number of sampled from a 2d distribution consistent with the input distribution. 
#COULD BE IMPROVED BY A RECTLINEAR BIVARIATE SPLINE
import sample2dDistribution #.sample2dDistribution import sample2dDistribution
key2 = list(mDFL.keys())[10]
key1 = list(mDFL.keys())[1]
EvsECCENsampler = sample2dDistribution.sample2dDistribution(counts=mDFL[key1][key2]['HdistNORM'],xcent=mDFL[key1][key2]['xcent'],ycent=mDFL[key1][key2]['ycent'])
#EvsECCENsampler.plotPMF()
#EvsECCENsampler.plotCMF()
nsamples = 10**4.
xind, yind, x, y = EvsECCENsampler.discreteInverseTransformSampler(nsamples)
EvsECCENsampler.plot2dPMF(x,y)

plt.close(65183)
fig = plt.figure(65183)
plt.scatter(x, y, alpha=0.1)
plt.xlabel('x cent')
plt.ylabel('y cent')
plt.show(block=False)

plt.draw()
####################################################################################################################

#### Function for marginalizing over joint probability distribution of y given x
def calc_fdY(Y, xmin, xmax, EVPOCpdf, xnew):
    """Calculates probability density of Y by integrating over X - as done in BrownCompleteness
    
    Args:
        Y (float ndarray):
            Value of Y
        xmin (float ndarray):
            Value of minimum projected separation (AU) from instrument
        xmax (float ndarray):
            Value of maximum projected separation (AU) from instrument
        EVPOCpdf () - pdf function as done in BrownCompleteness
        xnew () - as done in BrownCompleteness

    Returns:
        float:
            Value of probability density
    
    """
    
    f = np.zeros(len(xmin))
    for k, dY in enumerate(Y):
        f[k] = interpolate.InterpolatedUnivariateSpline(xnew,EVPOCpdf(xnew,dY),ext=1).integral(xmin[k],xmax[k])
        
    return f

#### Calculate the probability distribution of r ####################################
# FIRST plot a data cube for all values of r, since r(a,e,E)
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from matplotlib import ticker
eccen_rplot = np.linspace(start=ax2[0][1].get_xlim()[0],stop=ax2[0][1].get_xlim()[1],num=300)
E_rplot = np.linspace(start=ax2[0][10].get_xlim()[0],stop=ax2[0][10].get_xlim()[1],num=300)
SMA_rplot = np.linspace(start=ax2[0][6].get_xlim()[0],stop=ax2[0][6].get_xlim()[1],num=300)
eccen_x, SMA_x, E_x = np.meshgrid(eccen_rplot, SMA_rplot, E_rplot)

#Calculate all r
r = SMA_x*(1.-eccen_x**2.)/(1.+eccen_x*((np.cos(E_x)-eccen_x)/(1.-eccen_x*np.cos(E_x))))

plt.close(35483)
fig = plt.figure(35483,figsize=(6,10))
ax= fig.add_subplot(111, projection= '3d')
alpha = 0.5
surf1 = ax.contourf(eccen_rplot, SMA_rplot, r[:,:,0], zdir='z', offset=E_rplot[0],cmap='bwr',vmin=np.min(r), vmax=np.max(r),alpha=alpha)#, norm=mpl.colors.LogNorm())#,locator=ticker.LogLocator())#'afmhot')
surf15 = ax.contourf(eccen_rplot, SMA_rplot, r[:,:,50], zdir='z', offset=E_rplot[50],cmap='bwr',vmin=np.min(r), vmax=np.max(r),alpha=alpha)#, norm=mpl.colors.LogNorm())#,locator=ticker.LogLocator())#'afmhot')
surf2 = ax.contourf(eccen_rplot, SMA_rplot, r[:,:,100], zdir='z', offset=E_rplot[100],cmap='bwr',vmin=np.min(r), vmax=np.max(r),alpha=alpha)#, norm=mpl.colors.LogNorm())#,locator=ticker.LogLocator())#'afmhot')
surf25 = ax.contourf(eccen_rplot, SMA_rplot, r[:,:,150], zdir='z', offset=E_rplot[150],cmap='bwr',vmin=np.min(r), vmax=np.max(r),alpha=alpha)#, norm=mpl.colors.LogNorm())#,locator=ticker.LogLocator())#'afmhot')
surf3 = ax.contourf(eccen_rplot, SMA_rplot, r[:,:,200], zdir='z', offset=E_rplot[200],cmap='bwr',vmin=np.min(r), vmax=np.max(r),alpha=alpha)#, norm=mpl.colors.LogNorm())#,locator=ticker.LogLocator())#'afmhot')
surf35 = ax.contourf(eccen_rplot, SMA_rplot, r[:,:,250], zdir='z', offset=E_rplot[250],cmap='bwr',vmin=np.min(r), vmax=np.max(r),alpha=alpha)#, norm=mpl.colors.LogNorm())#,locator=ticker.LogLocator())#'afmhot')
surf4 = ax.contourf(eccen_rplot, SMA_rplot, r[:,:,-1], zdir='z', offset=E_rplot[-1],cmap='bwr',vmin=np.min(r), vmax=np.max(r),alpha=alpha)#, norm=mpl.colors.LogNorm())#,locator=ticker.LogLocator())#'afmhot')
ax.set_xlabel('e', weight='bold')
ax.set_ylabel('a in AU', weight='bold')
ax.set_zlabel('E in rad', weight='bold')
ax.set_xlim([0.,1.])
ax.set_ylim([np.min(SMA_rplot),np.max(SMA_rplot)])
ax.set_zlim([0.,2.*np.pi])
ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([0.75, 0.75, 1.2, 1])) # np.diag([0.5, 0.5, 1.0, 1]))
ax.elev = ax.elev*0.7

N=100
cmap = plt.get_cmap('bwr',N)
norm = mpl.colors.Normalize(vmin=np.min(r), vmax=np.max(r))
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, boundaries=np.arange(0.,np.ceil(np.max(r)),1.))
#DELETE cbar.set_label('Spacecraft ' + r'|\underline{r}|' + ' in AU', rotation=270)
cbar.set_label('Spacecraft r in AU', rotation=270, labelpad=30)
plt.tight_layout()
plt.show(block=False)
del eccen_rplot, SMA_rplot, E_rplot #these values should not be used elsewhere
######

#### Calculate probability densities for each point
#Calculate probability of e_i, a_i given e_i, E_i given e_i

#Define keys
ikey1 = 10 # Eccentric Anomaly Row
ikey2 = 1 #eccentricity column
key1 = list(pdData.keys())[ikey1] #ok #Eccentric\nAnomaly\n(rad)
key2 = list(pdData.keys())[ikey2] #ok #Eccentricity

#Calculate value ranges
eccen_range = np.linspace(start=np.min(pdData['Eccentricity']),stop=np.max(pdData['Eccentricity']),num=300, endpoint=True)
E_range = np.linspace(start=np.min(pdData['Eccentric\nAnomaly\n(rad)']),stop=np.max(pdData['Eccentric\nAnomaly\n(rad)']),num=300, endpoint=True)

#Verify the intgral over the eccentricity gaussian KDE is 1.
eccen_normalization_constant = kdes[ikey2].integrate_box_1d(-np.inf, np.inf)

#Calculate probabilities for each eccen in the r data cube
P_eccen_rplot = kdes[ikey2](eccen_range)

#Demonstrates Hdist sums to total # of SC
print(str(np.sum(mDFL[key1][key2]['Hdist'])) + ' : this value should be approximately 17966 spacecraft')
#Demonstrates HdistNORM does approximately sum to 1
xaxismin, xaxismax, yaxismin, yaxismax = get_pdDataLimits(pdData,key2,key1)
print(str(np.sum(mDFL[key1][key2]['HdistNORM'])*(xaxismax - xaxismin)*(yaxismax - yaxismin)) + ' : this value should be approximately Number of SC 17966')
#DELETE print(str(np.sum(mDFL[key1][key2]['HdistNORM'])*(np.max(mDFL[key1][key2]['xcent']) - np.min(mDFL[key1][key2]['xcent']))*(np.max(mDFL[key1][key2]['ycent']) - np.min(mDFL[key1][key2]['ycent']))) + ' : this value should be approximately Number of SC 17966')
#Demonstrates EVPOCpdf integral sums to 1
EVPOCpdf = mDFL[key1][key2]['EVPOCpdf']
xaxismin, xaxismax, yaxismin, yaxismax = get_pdDataLimits(pdData,key2,key1)
print(str(EVPOCpdf.integral(xaxismin,xaxismax,yaxismin,yaxismax)) + ' : this value should be approximately 1')
print(str(EVPOCpdf.integral(-np.inf,np.inf,-np.inf,np.inf)) + ' : this value should be approximately 1')


num=len(mDFL[key1][key2]['ycent'])#1000
Es = mDFL[key1][key2]['ycent']
#eccen_center = sampler[ikey2](1) + np.zeros(num) # randomly select eccentricity center value
eccen_center = 0.0005+ np.zeros(num) # randomly select eccentricity center value
deccen = (np.max(pdData['Eccentricity']) - np.min(pdData['Eccentricity']))/1000. # this is arbitrarily chosen, could be smaller 

#This figure demonstrates the edge effects of the RBS and fdY calculation process 
plt.close(6565984397877)
fign = plt.figure(6565984397877)
E_marginalized_eccen = calc_fdY( Es, eccen_center-deccen/2., eccen_center+deccen/2., mDFL[key1][key2]['EVPOCpdf'], mDFL[key1][key2]['xnew'])
plt.plot(Es,E_marginalized_eccen,color='blue')
Es = np.linspace(start=np.min(pdData['Eccentric\nAnomaly\n(rad)']),stop=np.max(pdData['Eccentric\nAnomaly\n(rad)']),num=num,endpoint=True)
E_marginalized_eccen = calc_fdY( Es, eccen_center-deccen/2., eccen_center+deccen/2., mDFL[key1][key2]['EVPOCpdf'], mDFL[key1][key2]['xnew'])
plt.plot(Es,E_marginalized_eccen,color='red')
plt.xlabel('E|e_i')
plt.ylabel('Probability(E|e_i)')
plt.show(block=False)


#E_marginalized_eccen = mDFL[key1][key2]['EVPOCpdf'].integral()


###########################################################################
#Calculates Probability of E given eccen
eccen_grid, E_grid = np.meshgrid(eccen_range,E_range) 
P_E_given_eccen = mDFL[key1][key2]['vectorizedEVPOCpdf'](eccen_grid,E_grid)
#NOTE: WHILE INTERESTING, THIS IS NOT USEFUL BECAUSE A SINGLE ECCENTRICITY IS GENERATED SO ALL ECCENTRICITIES MUST BE THE SAME
#Calculate P(E|e)*P(e)
for i in np.arange(len(eccen_range)):
    P_E_given_eccen[i,:] = P_eccen_rplot[i]*P_E_given_eccen[i,:]
#Normalize P(E|e)*P(e)
P_E_given_eccen = P_E_given_eccen/np.sum(P_E_given_eccen)
plt.close(321683687)
fig = plt.figure(321683687,figsize=(8,6))
#ax= fig.add_subplot(111)#, projection= '3d')
alpha = 0.5
plt.contourf(eccen_grid, E_grid, P_E_given_eccen)
plt.xlabel('e', weight='bold')
plt.ylabel('E in rad', weight='bold')
#ax.set_zlabel('Probability')
plt.colorbar(cmap='bwr')
plt.show(block=False)
#####################################################################################



# #### Calculate Rectlinear bivariate splines for the data #######################################################
# from scipy import interpolate

# #Just for Eccen vs SMA
# xmin = ax2[1][6].get_xlim()[0]
# xmax = ax2[1][6].get_xlim()[1]
# xedges = np.linspace(start=xmin,stop=xmax,num=100)
# xcent = 0.5*(xedges[1:]+xedges[:-1])
# ymin = ax2[1][6].get_ylim()[0]
# ymax = ax2[1][6].get_ylim()[1]
# yedges = np.linspace(start=ymin,stop=ymax,num=100)
# ycent = 0.5*(yedges[1:]+yedges[:-1])
# xnew = np.hstack((0.0,xcent,xmax))
# ynew = np.hstack((ymin,ycent,ymax))

# Cpdf = np.pad(Cpdf,1,mode='constant')

# #save interpolant to object
# self.Cpdf = Cpdf
# self.EVPOCpdf = interpolate.RectBivariateSpline(xnew, ynew, Cpdf.T)



#### Calculate Visual Magnitude of Spacecraft ######################################################################
####################################################################################################################





