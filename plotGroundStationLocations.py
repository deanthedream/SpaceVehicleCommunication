""" Plot Ground Station Locations
Written By: Dean Keithly
Written On: 6/11/2019
"""

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import os
import csv

#### Load in Ground Station Locations from Vallado 2011 "Simulating Space Survelliance Networks"
GSpaperFilePath = '/home/dean/Documents/AFRL2019'
fileNames = ['GroundStationDataTable1.csv', # US AF Space Surveillance Network Sensors
'GroundStationDataTable2.csv', # US AF Space Surveillance Network Sensors
'GroundStationDataTable3.csv', # Russian Space Surveillance Sensors
'GroundStationDataTable4.csv', # European Space Surveillance Sensors
'GroundStationDataTable5.csv', # United States Air Force Satellite Control Network
'GroundStationDataTable6.csv', #Chinese Space Surveillance System
'GroundStationDataTable7.csv', #Other Governmental Space Surveillance Systems
'GroundStationDataTable8.csv', # International Scientific Optical Network Sensor
'GroundStationDataTable9.csv'] # Satellite Laser Ranging Sites
GSpaperData = list()
for i in np.arange(len(fileNames)):
    with open(os.path.join(GSpaperFilePath,fileNames[i]), 'r') as f:
        lineReader = csv.reader(f,delimiter=',')
        tmp = list()
        for line in lineReader:
            tmp.append(line)
        GSpaperData.append(tmp)
"""
GSpaperData[fileNameIndex][len(number of stations)][columns]
columns are:
    0:"ID #" the station ID (sometimes multiplied by 10)
    1:"Location" the name of the city and/or country
    2:"Name"  Name of the place
    3:"Type" type of observation station (typically Opt)
    4:"Latitude" lat in deg
    5:"Longitude" lon in deg
    6:"Alt (m)" altitude in m above sea-level
    7:"Open" date opened
    8:"Close" date closed (mostly empty)
    9:"Notes" uncoordinated notes
In cases where the ID # were unknown, radars were given “555” prefixes,
phased arrays a “666” prefix, and optical sensors were given a
“999” prefix within each sensor network.
For all other ID # in the thousands, the ID #'s' were multiplied by 10
if they represent the central location of a sensor cluster
"""


colors =  ['b','b','r','g','b','r','grey',  'purple',   'orange']
markers = ['o','o','p','d','d','d','o',     's',        'x']


# lon_0 is central longitude of projection.
# resolution = 'c' means use crude resolution coastlines.
m = Basemap(projection='hammer',lon_0=0,resolution='c')
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))

for i in np.arange(len(GSpaperData)):
    for j in np.arange(len(GSpaperData[i])):
        if j == 0:
            continue
        lat = float(GSpaperData[i][j][4])
        lon = float(GSpaperData[i][j][5])
        x,y = m(lon, lat)
        m.plot(x, y, color=colors[i], marker=markers[i], markersize=5)

m.drawmapboundary(fill_color='aqua')
plt.title("Hammer Projection")
plt.show(block=False)

