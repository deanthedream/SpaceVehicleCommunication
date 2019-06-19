#UDLsensorGroundStations.py

#This replicates UDL querying from UDLqueries to collect
#a list of sensors (that I have access to) and extracts their lon/lat where available

import numpy as np
from grabUDL import grabUDL as grabUDL
import astropy.units as u
import getDateTime

g = grabUDL()


###### Query All Sensors #####################################################
url_allSensors = "https://unifieddatalibrary.com/udl/sensor"
dat_allSensors, col_allSensors = g.queryUDL(url_allSensors)
arr_allSensors = dat_allSensors[col_allSensors]#.values # extracts data as an array

for key in ['classificationMarking', 'dataControl', 'dataMode', 'entity',
       'idSensor', 'sensorName', 'sensorNumber', 'sensorObservationType',
       'sensorType', 'sensorcharacteristics', 'sensorlimitsCollection',
       'shortName', 'source', 'taskable']:
    print('key: ' + key + ' data: ' + str(arr_allSensors[key][736]))

##### Extracts list of lon-lat locations of all queried Sensors
lon = list()
lat = list()
for i in np.arange(len(arr_allSensors['entity'])):
    #All sensors have location field
    #print(str(i) + ': ' + str(arr_allSensors['entity'][i]['location']))
    if 'lat' in arr_allSensors['entity'][i]['location']:
        lat.append(arr_allSensors['entity'][i]['location']['lat'])
        lon.append(arr_allSensors['entity'][i]['location']['lon'])
        print(str(i) + ', lat=' + str(arr_allSensors['entity'][i]['location']['lat']) + 
            ', lon=' + str(arr_allSensors['entity'][i]['location']['lon']))
UDLsensorData = {'lon':lon, 'lat':lat}
#########################################################################################

def plotUDLsensorsGS_basemap(m,UDLsensorData):
    """
    """
    for i in np.arange(len(UDLsensorData['lon'])):
        lat = float(UDLsensorData['lat'][i])
        lon = float(UDLsensorData['lon'][i])
        x,y = m(lon, lat)
        m.plot(x, y, color='k', marker='o', markersize=5, alpha=0.5)
        #lat = float(UDLsensorData[i][1][4])
        #lon = float(UDLsensorData[i][1][5])
        #x,y = m(lon, lat)
        #m.plot(x, y, color=color[i], linestyle=None, marker=marker[i], markersize=5, label=networkLabels[i])
    return m


