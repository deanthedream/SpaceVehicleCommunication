# UDL queries

import numpy as np
from grabUDL import grabUDL as grabUDL
import astropy.units as u

g = grabUDL()


#Example Query
url_example = "https://unifieddatalibrary.com/udl/elset?epoch=%3E2018-11-01T00:00:00.000000Z&satNo=25544"
dat_example, col_example = g.queryUDL(url_example)

#Query All Sensors
url_allSensors = "https://unifieddatalibrary.com/udl/sensor"
dat_allSensors, col_allSensors = g.queryUDL(url_allSensors)
arr_allSensors = dat_allSensors[col_allSensors]#.values # extracts data as an array


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




### Search for logical earth position indicators in sensors data
possibleLocNames0 = [param['name'] for param in dat_senorshelp['parameters'] if param['unitOfMeasure'] == 'degrees']
paramFields = [param for param in dat_senorshelp['parameters'] if param['name'] in ['lat',
 'lon', 'leftClockAngle', 'rightClockAngle', 'boresight', 'boresightOffAngle', 'maxDeviationAngle', 'leftGeoBeltLimit', 'rightGeoBeltLimit']]

sensorLats = dat_allSensors['lat']
sensorLons = dat_allSensors['lon']

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