""" Plot Ground Station Locations
Written By: Dean Keithly
Written On: 6/11/2019
"""

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import vallado2014GroundStations # loads data extracted from Vallado2014
import UDLsensorGroundStations # contains a query to all sensors to ellicit ground station locations

# lon_0 is central longitude of projection.
# resolution = 'c' means use crude resolution coastlines.
plt.figure(1,figsize=(10,6))
m = Basemap(projection='hammer',lon_0=0,resolution='c')
m.drawcoastlines()
#m.fillcontinents(color='coral',lake_color='aqua')
m.fillcontinents(color='white',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))

#Plot Vallado Ground Station Locations
m = vallado2014GroundStations.plotVallado2014GS_basemap(m,vallado2014GroundStations.GSpaperData)

#Plot Sensor Ground Station Locations
m = UDLsensorGroundStations.plotUDLsensorsGS_basemap(m,UDLsensorGroundStations.UDLsensorData)

m.drawmapboundary(fill_color='aqua')
plt.legend(loc='lower center', ncol=5, bbox_to_anchor=(0.5, -0.1),
    fancybox=True, shadow=True)
plt.title("Vallado Ground Station Locations")
plt.tight_layout()
plt.show(block=False)




#POTENTIAL FUTURE IMPROVEMENT
"""
I had luck installing Cesium - looks a lot like the STK - AGI Earth viewing tool
Be sure to run cd ~/Cesium and node server.js
I could load the python module and do the following example from https://github.com/sinhrks/cesiumpy
import cesiumpy
v=cesiumpy.Viewer()
v.entities.add(cesiumpy.Box(dimensions=(40e4, 30e4, 50e4),material=cesiumpy.color.RED, position=(-120, 40, 0)))
v #In a jupyter notebook, this would display something (supposedly)
v.to_html() #outputs <script> stuff to be put into an HTML but I can't get this to work right...

Here are some links
https://cesium.com/docs/tutorials/getting-started/
https://cesiumpy.readthedocs.io/en/latest/
https://cesiumpy.readthedocs.io/en/latest/basics.html
Note I created an account to get a 'token' which can be found at https://cesium.com/ion/tokens

An interesting example scraping locations from wikipedia
https://cesiumpy.readthedocs.io/en/latest/examples.html

#dont know why but here is some pysat stuff https://buildmedia.readthedocs.org/media/pdf/pysat/master/pysat.pdf
"""

