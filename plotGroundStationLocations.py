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
fileNames = ['GroundStationDataTable1.csv','GroundStationDataTable2.csv','GroundStationDataTable3.csv',\
'GroundStationDataTable4.csv','GroundStationDataTable5.csv','GroundStationDataTable6.csv',\
'GroundStationDataTable7.csv','GroundStationDataTable8.csv','GroundStationDataTable9.csv']
GSpaperData = list()
for i in np.arange(len(fileNames)):
    with open(os.path.join(GSpaperFilePath,fileNames[i]), 'r') as f:
        GSpaperData.append(csv.reader(f,delimiter=','))







# lon_0 is central longitude of projection.
# resolution = 'c' means use crude resolution coastlines.
m = Basemap(projection='hammer',lon_0=0,resolution='c')
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
m.drawmapboundary(fill_color='aqua')
plt.title("Hammer Projection")
plt.show(block=False)

