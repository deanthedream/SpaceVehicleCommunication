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
########################################################################

#### Plot Scatter Matrix ###############################################
from pandas.plotting import scatter_matrix
import pandas as pd 
#DELETE from matplotlib.ticker import StrMethodFormatter

#Format data to turn into pandas frame
pdData = {'Mean\nMotion\n(rad/min)':no, 'Eccentricity':ecco, 'Inclination\n(rad)':inclo,\
    'Mean\nAnomaly\n(rad)':mo, 'Arg. of\nPerigee\n(rad)':argpo, 'Longitude\nAscending\nNode\n(rad)':nodeo,\
    'Semi-major\nAxis\n(Earth Radii)':a,'Apoapsis\nAltitude\n(Earth Radii)':alta, 'Periapsis\nAtlitude\n(Earth radii)':altp, 'Time Since\nEpoch (JD)':jdsatepoch}

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
ax2 = scatter_matrix(df, alpha=0.05, diagonal='kde', ax=ax, **{'color':'black'})#, **kwds)
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
                ax_sub2.set_yticklabels(tmplabels)
        label2 = ax_sub2.get_xlabel()
        ax_sub2.set_xlabel(label2, weight='bold')
plt.show(block=False)
print('Done Plotting Scatter Matrix')

#OK### Plot Eccen. vs SMA relation based on R_earth Min
#EXCLUSION REGION IS TOP LEFT
SMA = np.linspace(start=ax2[1][6].get_xlim()[0],stop=ax2[1][6].get_xlim()[1],num=30) #Semi-Major Axis
R_earth = 1. #In earth radii
minAlt = 0.#LEO min alt is technically earth surface 400./6356. #assumes 400km is min LEO alt, Earth Radius is 6356 km from Wikipedia
R = R_earth + minAlt
ECCEN = 1. - R/SMA
ax2[1][6].plot(SMA,ECCEN,color='b',alpha=0.4)
plt.show(block=False)
#OK### Plot SMA vs Eccen. relation based on R_earth Min
#EXCLUSION REGION IS BOTTOM RIGHT
ECCEN = np.linspace(start=ax2[6][1].get_xlim()[0],stop=ax2[6][1].get_xlim()[1],num=30) #Semi-Major Axis
R_earth = 1. #In earth radii
minAlt = 0.#LEO min alt is technically earth surface 400./6356. #assumes 400km is min LEO alt, Earth Radius is 6356 km from Wikipedia
R = R_earth + minAlt
SMA = R/(1.-ECCEN)
ax2[6][1].plot(ECCEN,SMA,color='b',alpha=0.4)
plt.show(block=False)
#OK### Plot Eccen. vs Apoapsis relation based on R_earth Min
#EXCLUSION REGION IS TOP LEFT
APOAPSIS = np.linspace(start=ax2[1][7].get_xlim()[0],stop=ax2[1][7].get_xlim()[1],num=30) #Semi-Major Axis
R_earth = 1. #In earth radii
minAlt = 0.#LEO min alt is technically earth surface 400./6356. #assumes 400km is min LEO alt, Earth Radius is 6356 km from Wikipedia
R = R_earth + minAlt
ECCEN = 2.*APOAPSIS/(R+APOAPSIS) - 1.
ax2[1][7].plot(APOAPSIS-R_earth,ECCEN,color='b',alpha=0.4)
plt.show(block=False)
#OK### Plot Apoapsis vs Eccen. relation based on R_earth Min
#EXCLUSION REGION IS BOTTOM RIGHT
ECCEN = np.linspace(start=ax2[7][1].get_xlim()[0],stop=ax2[7][1].get_xlim()[1],num=30) #eccentricity
R_earth = 1. #In earth radii
minAlt = 0.#LEO min alt is technically earth surface 400./6356. #assumes 400km is min LEO alt, Earth Radius is 6356 km from Wikipedia
R = R_earth + minAlt
APOAPSIS = R*(ECCEN+1.)/(1.-ECCEN)
ax2[7][1].plot(ECCEN,APOAPSIS-R_earth,color='b',alpha=0.4)
plt.show(block=False)
#OK### Plot APOAPSIS vs SMA
SMA = np.linspace(start=ax2[7][6].get_xlim()[0],stop=ax2[7][6].get_xlim()[1],num=30) #Semi-Major Axis
APOAPSIS = SMA-R
ax2[7][6].plot(SMA,APOAPSIS,color='r',alpha=0.4)
plt.show(block=False)
#OK### Plot SMA vs APOAPSIS
APOAPSIS = np.linspace(start=ax2[6][7].get_xlim()[0],stop=ax2[6][7].get_xlim()[1],num=30) #Semi-Major Axis
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
ECCEN = np.linspace(start=ax2[8][1].get_xlim()[0],stop=ax2[8][1].get_xlim()[1],num=30) #ECCEN RANGE
ax2[8][1].plot(ECCEN,R-R + np.zeros(len(ECCEN)),color='r',alpha=0.4)
plt.show(block=False)
#OK### Plot Minimum APOAPSIS vs PERIAPSIS, since APOAPSIS must be larger than PERIAPSIS #BOTTOM RIGHT EXCLUSION ZONE
PERIAPSIS = np.linspace(start=ax2[7][8].get_xlim()[0],stop=ax2[7][8].get_xlim()[1],num=30) #PERIAPSIS raNge
APOAPSIS = PERIAPSIS
ax2[7][8].plot(PERIAPSIS,APOAPSIS,color='r',alpha=0.4)
plt.show(block=False)
#OK### Plot MAXIMUM PERIAPSIS vs APOAPSIS, since APOAPSIS must be larger than PERIAPSIS #TOP LEFT EXCLUSION ZONE
APOAPSIS = np.linspace(start=ax2[8][7].get_xlim()[0],stop=ax2[8][7].get_xlim()[1],num=30) #APOAPSIS raNge
PERIAPSIS = APOAPSIS
ax2[8][7].plot(APOAPSIS,PERIAPSIS,color='r',alpha=0.4)
plt.show(block=False)
#OK### Plot Mean Motion vs SMA
SMA = np.linspace(start=ax2[0][6].get_xlim()[0],stop=ax2[0][6].get_xlim()[1],num=30) #Semi-Major Axis
mu = const.GM_earth
MEANMOTION = np.sqrt(mu.value*(60.*60.)/(const.R_earth.value**3.)/SMA**3.)
ax2[0][6].plot(SMA,MEANMOTION,color='r',alpha=0.4)
plt.show(block=False)
#OK### Plot SMA vs Mean Motion
MEANMOTION = np.linspace(start=ax2[6][0].get_xlim()[0],stop=ax2[6][0].get_xlim()[1],num=30) #Mean Motion
mu = const.GM_earth
MEANMOTION = np.sqrt(mu.value*(60.*60.)/(const.R_earth.value**3.)/SMA**3.)
SMA = (mu.value*(60.*60.)/(const.R_earth.value**3.)/MEANMOTION**2.)**(1./3.)
ax2[6][0].plot(MEANMOTION,SMA,color='r',alpha=0.4)
plt.show(block=False)
#OK### Plot Mean Motion vs APOAPSIS
R_earth = 1.
APOAPSIS = np.linspace(start=ax2[0][7].get_xlim()[0],stop=ax2[0][7].get_xlim()[1],num=30) #APOAPSIS
MEANMOTION = np.sqrt(mu.value*(60.*60.)/(const.R_earth.value**3.)/(APOAPSIS)**3.)
ax2[0][7].plot(APOAPSIS-R_earth,MEANMOTION,color='r',alpha=0.4) #should be a lower bound
MEANMOTION = np.sqrt(8.*mu.value*(60.*60.)/(const.R_earth.value**3.)/(APOAPSIS)**3.) #should be an upper limit
ax2[0][7].plot(APOAPSIS-R_earth,MEANMOTION,color='r',alpha=0.4)
plt.show(block=False)
#NOT TOTALLY SURE ABOUT THE X-AXIS VALUES. NEED TO DOUBLE CHECK WHETHER APOAPSIS HERE IS CORRECT
#OK### Plot Mean Motion vs PERIAPSIS
R_earth = 1.
PERIAPSIS = np.linspace(start=ax2[0][8].get_xlim()[0],stop=ax2[0][8].get_xlim()[1],num=30) #PERIAPSIS
ax2[0][8].plot(PERIAPSIS-R_earth,0. + np.zeros(len(PERIAPSIS)),color='r',alpha=0.4) # since 1-e_max = 0 
MEANMOTION = np.sqrt(mu.value*(60.*60.)/(const.R_earth.value**3.)/(PERIAPSIS)**3.)
ax2[0][8].plot(PERIAPSIS-R_earth,MEANMOTION,color='r',alpha=0.4) # since 1-e_min = 0 
plt.show(block=False)



#### Save Scatter Data
import datetime
import re
import os
# Save to a File
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
###########################################################################

#comment-out to improve speed
#plotOrbitDist_scatterMatrix(pdData)
#############################################################################

#### Plot Spacecraft Parameters in "contour matrix" #######################
from plotContourFromScatter import plotContourFromScatter
from plotContourFromScatter import plotKDEfromScatter

mDFL = {} # Create Struct to hold all plot objects
for key1 in pdData.keys():
    mDFL[key1] = {}

#Create all axis options
fignum=516841
plt.close(fignum)
plt.rc('axes',linewidth=2)
plt.rc('lines',linewidth=2)
plt.rcParams['axes.linewidth']=2
plt.rc('font',weight='bold')
outAX = plt.subplots(10,10, num=fignum, figsize=(12,12))
outAX[0].subplots_adjust(bottom=0.15)
#outAX[0] is the figure object
for ikey1 in np.arange(len(pdData.keys())):
    for ikey2 in np.arange(len(pdData.keys())):
        key1 = list(pdData.keys())[ikey1]
        key2 = list(pdData.keys())[ikey2]
        if key1 == key2:
            outAX[1][ikey1][ikey2] = plotKDEfromScatter(pdData[key1], ax=outAX[1][ikey1][ikey1])
            #mDFL[key1][key2] = plt.gca()
        else:
            outAX[1][ikey1][ikey2] = plotContourFromScatter(pdData[key1],pdData[key2],ax=outAX[1][ikey1][ikey2],nbins=11,figsize=(6,8))
            #mDFL[key1][key2] = plt.gca()
        if not key2 == list(pdData.keys())[0]: # the y-axis keys
            outAX[1][ikey1][ikey2].get_yaxis().set_visible(False)
        else:
            outAX[1][ikey1][ikey2].set_ylabel('\n' + list(pdData.keys())[ikey1],weight='bold',rotation=0,labelpad=40) #set ylabel

            #NEED ADJUSTMENT IF KEY2=0, USE Y AXIS LABELS OF PLOT TO RIGHT
        if not key1 == list(pdData.keys())[-1]: # the x-axis keys
            outAX[1][ikey1][ikey2].get_xaxis().set_visible(False)
        else:
            outAX[1][ikey1][ikey2].tick_params(axis='x', rotation=90)
            outAX[1][ikey1][ikey2].set_xlabel(list(pdData.keys())[ikey2],weight='bold') #set xlabel
            
plt.subplots_adjust(wspace=0, hspace=0)
plt.show(block=False)
print('done creating all contour plots')
############################################################################################################


#### Calculate PDF and InverseTransform Sampler From data #################################################################################
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

    xmin = ax2[i][i].get_xlim()[0]
    xmax = ax2[i][i].get_xlim()[1]
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


#### Calculate Visual Magnitude of Spacecraft ######################################################################

####################################################################################################################





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


