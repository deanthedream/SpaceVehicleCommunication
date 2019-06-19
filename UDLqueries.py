# UDL queries

import numpy as np
from grabUDL import grabUDL as grabUDL
import astropy.units as u
import getDateTime

g = grabUDL()


#Example Query
url_example = "https://unifieddatalibrary.com/udl/elset?epoch=%3E2018-11-01T00:00:00.000000Z&satNo=25544"
dat_example, col_example = g.queryUDL(url_example)

#### QUERY HELP COMPONENTS ##################################
#Query Help All Sensors
url_senorshelp = "https://unifieddatalibrary.com/udl/sensor/queryhelp"
dat_senorshelp, col_senorshelp = g.queryUDL(url_senorshelp)
dat_senorshelp['parameters'] = g.addAstropyUnits_to_queryHelp(dat_senorshelp)

#query observation help
url_observationhelp = "https://unifieddatalibrary.com/udl/observation/queryhelp"
dat_observationhelp, col_observationhelp = g.queryUDL(url_observationhelp)
dat_observationhelp['parameters'] = g.addAstropyUnits_to_queryHelp(dat_observationhelp)

#query rfobservation help
url_rfobservationhelp = "https://unifieddatalibrary.com/udl/rfobservation/queryhelp"
dat_rfobservationhelp, col_rfobservationhelp = g.queryUDL(url_rfobservationhelp)
dat_rfobservationhelp['parameters'] = g.addAstropyUnits_to_queryHelp(dat_rfobservationhelp)

#query eoobservation help
url_eoobservationhelp = "https://unifieddatalibrary.com/udl/eoobservation/queryhelp"
dat_eoobservationhelp, col_eoobservationhelp = g.queryUDL(url_eoobservationhelp)
dat_eoobservationhelp['parameters'] = g.addAstropyUnits_to_queryHelp(dat_eoobservationhelp)

#query radarobservation help
url_radarobservationhelp = "https://unifieddatalibrary.com/udl/radarobservation/queryhelp"
dat_radarobservationhelp, col_radarobservationhelp = g.queryUDL(url_radarobservationhelp)
dat_radarobservationhelp['parameters'] = g.addAstropyUnits_to_queryHelp(dat_radarobservationhelp)

#query elset help - on orbit parameter querying
url_elsethelp = "https://unifieddatalibrary.com/udl/elset/queryhelp"
dat_elsethelp, col_elsethelp = g.queryUDL(url_elsethelp)
dat_elsethelp['parameters'] = g.addAstropyUnits_to_queryHelp(dat_elsethelp)
##############################################################################

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
        lon.append(arr_allSensors['entity'][i]['location']['lat'])
        lat.append(arr_allSensors['entity'][i]['location']['lon'])
        print(str(i) + ', lat=' + str(arr_allSensors['entity'][i]['location']['lat']) + 
            ', lon=' + str(arr_allSensors['entity'][i]['location']['lon']))


#Extract set of sensor IDs
#Valid keys to dat_allSensors
# ['classificationMarking', 'dataControl', 'dataMode', 'entity',
#        'idSensor', 'sensorName', 'sensorNumber', 'sensorObservationType',
#        'sensorType', 'sensorcharacteristics', 'sensorlimitsCollection',
#        'shortName', 'source', 'taskable']
set_sensorIDS = arr_allSensors['idSensor']
assert len(arr_allSensors['idSensor']) == len(set(arr_allSensors['idSensor'])), 'Error: There is a duplicate sensor'
#print(len(set(arr_allSensors['sensorNumber']))) #this indicates not all sensor numbers are unique. There are multiple nan, but there are also multiple nan in the set... interesting
##############################################################################

#### Test Query specific sensor ############################################
#ADS_RiverlandDingo_03
#url_ADS_RiverlandDingo_03 = "https://unifieddatalibrary.com/udl/sensor/ADS_RiverlandDingo_03"
#url_ADS_RiverlandDingo_03 = g.addQueryParameter(url_ADS_RiverlandDingo_03, 'idSensor', '=', 'ADS_RiverlandDingo_03')
#dat_ADS_RiverlandDingo_03, col_ADS_RiverlandDingo_03 = g.queryUDL(url_ADS_RiverlandDingo_03)
#t1, t2 = g.queryUDL("https://unifieddatalibrary.com/udl/sensor/ADS_RiverlandDingo")

#ADS_RiverlandDingo_03
url_ADS_RiverlandDingo_03 = "https://unifieddatalibrary.com/udl/observation"
#url_ADS_RiverlandDingo_03 = g.addQueryParameter(url_ADS_RiverlandDingo_03, 'id', '=', '5')
#url_ADS_RiverlandDingo_03 = g.addQueryParameter(url_ADS_RiverlandDingo_03, 'idSensor', '=', 'ADS_RiverlandDingo_03')
url_ADS_RiverlandDingo_03 = g.addQueryParameter(url_ADS_RiverlandDingo_03, 'obTime', '=>', getDateTime.formUDLtimeString(getDateTime.getULD_UTC(-200.0)))
url_ADS_RiverlandDingo_03 = g.addQueryParameter(url_ADS_RiverlandDingo_03, 'maxResults', '=', 3)
print(url_ADS_RiverlandDingo_03)
dat_ADS_RiverlandDingo_03, col_ADS_RiverlandDingo_03 = g.queryUDL(url_ADS_RiverlandDingo_03)

print(saltyburrito)


url_1 = "https://unifieddatalibrary.com/udl/rfobservation"
#url_ADS_RiverlandDingo_03 = g.addQueryParameter(url_ADS_RiverlandDingo_03, 'id', '=', '5')
#url_ADS_RiverlandDingo_03 = g.addQueryParameter(url_ADS_RiverlandDingo_03, 'idSensor', '=', 'ADS_RiverlandDingo_03')
url_1 = g.addQueryParameter(url_1, 'obTime', '=>', getDateTime.formUDLtimeString(getDateTime.getULD_UTC(-40.0)))
url_1 = g.addQueryParameter(url_1, 'maxResults', '=', 3)
print(url_1)
dat_1, col_1 = g.queryUDL(url_1)
print('dat_1')
print(dat_1)

url_2 = "https://unifieddatalibrary.com/udl/radarobservation"
#url_ADS_RiverlandDingo_03 = g.addQueryParameter(url_ADS_RiverlandDingo_03, 'id', '=', '5')
#url_ADS_RiverlandDingo_03 = g.addQueryParameter(url_ADS_RiverlandDingo_03, 'idSensor', '=', 'ADS_RiverlandDingo_03')
url_2 = g.addQueryParameter(url_2, 'obTime', '=>', getDateTime.formUDLtimeString(getDateTime.getULD_UTC(-40.0)))
url_2 = g.addQueryParameter(url_2, 'maxResults', '=', 3)
print(url_2)
dat_2, col_2 = g.queryUDL(url_2)
print('dat_2')
print(dat_2)

url_3 = "https://unifieddatalibrary.com/udl/eoobservation"
#url_ADS_RiverlandDingo_03 = g.addQueryParameter(url_ADS_RiverlandDingo_03, 'id', '=', '5')
#url_ADS_RiverlandDingo_03 = g.addQueryParameter(url_ADS_RiverlandDingo_03, 'idSensor', '=', 'ADS_RiverlandDingo_03')
url_3 = g.addQueryParameter(url_3, 'obTime', '=>', getDateTime.formUDLtimeString(getDateTime.getULD_UTC(-40.0)))
url_3 = g.addQueryParameter(url_3, 'maxResults', '=', 3)
print(url_3)
dat_3, col_3 = g.queryUDL(url_3)
print('dat_3')
print(dat_3)

url_4 = "https://unifieddatalibrary.com/udl/sensor"
#url_ADS_RiverlandDingo_03 = g.addQueryParameter(url_ADS_RiverlandDingo_03, 'id', '=', '5')
#url_ADS_RiverlandDingo_03 = g.addQueryParameter(url_ADS_RiverlandDingo_03, 'idSensor', '=', 'ADS_RiverlandDingo_03')
url_4 = g.addQueryParameter(url_4, 'obTime', '=>', getDateTime.formUDLtimeString(getDateTime.getULD_UTC(-40.0)))
url_4 = g.addQueryParameter(url_4, 'maxResults', '=', 3)
print(url_4)
dat_4, col_4 = g.queryUDL(url_4)
print('dat_4')
print(dat_4)




print(saltyburrito)
#arr_aADS_RiverlandDingo_03 = dat_ADS_RiverlandDingo_03[col_ADS_RiverlandDingo_03]#.values # extracts data as an array

#####################################################



### Search for logical earth position indicators in sensors data
possibleLocNames0 = [param['name'] for param in dat_senorshelp['parameters'] if param['unitOfMeasure'] == 'degrees']
paramFields = [param for param in dat_senorshelp['parameters'] if param['name'] in ['lat',
 'lon', 'leftClockAngle', 'rightClockAngle', 'boresight', 'boresightOffAngle', 'maxDeviationAngle', 'leftGeoBeltLimit', 'rightGeoBeltLimit']]


#sensorLats = dat_allSensors['lat'] #doesn;t work
#sensorLons = dat_allSensors['lon']

out = [dat for dat in dat_allSensors['sensorcharacteristics'] if not dat == []]
#result unsuccessful. did not find lat and lon coords of sensors


### Search for logical earth position indicators in observation data
possibleLocNames1 = [param['name'] for param in dat_observationhelp['parameters'] if param['unitOfMeasure'] == 'degrees']
paramFields = [param for param in dat_observationhelp['parameters'] if param['name'] in ['azimuth', 'elevation', 'ra', 'declination']]


possibleLocNames2 = [param['name'] for param in dat_rfobservationhelp['parameters'] if param['unitOfMeasure'] == 'degrees']
#Query All rfobservation
url_allrfobservation = "https://unifieddatalibrary.com/udl/rfobservation"
dat_allrfobservation, col_allrfobservation = g.queryUDL(url_allrfobservation)
#arr_allrfobservation = dat_allrfobservation[col_allrfobservation]#.values # extracts data as an array
rfobservationLats2 = dat_allrfobservation['senlat']
rfobservationLons2 = dat_allrfobservation['senlon']


#NOPE
possibleLocNames3 = [param['name'] for param in dat_radarobservationhelp['parameters'] if param['unitOfMeasure'] == 'degrees']

#YEP
possibleLocNames4 = [param['name'] for param in dat_eoobservationhelp['parameters'] if param['unitOfMeasure'] == 'degrees']
#Query All rfobservation
url_alleoobservation = "https://unifieddatalibrary.com/udl/eoobservation"
dat_alleoobservation, col_alleoobservation = g.queryUDL(url_alleoobservation)
arr_alleoobservation = dat_alleoobservation[col_alleoobservation]#.values # extracts data as an array
eoobservationLats4 = dat_alleoobservation['senlat']
eoobservationLons4 = dat_alleoobservation['senlon']




#EXAMPLE QUERY
#https://unifieddatalibrary.com/udl/observation?obTime=%3E2019-06-01T00%3A00%3A00.000000Z