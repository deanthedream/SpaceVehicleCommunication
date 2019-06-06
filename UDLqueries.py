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


#Query Help All Sensors
url_help = "https://unifieddatalibrary.com/udl/sensor/queryhelp"
dat_help, col_help = g.queryUDL(url_help)

params = dat_help['parameters']

#Trying to convert these into manipulatable quantities in astropy
#{'units of magnitude', 'Hz', None, 'degrees', 'km/sec', 'rad/min', 'Watts', 'decibels', 'kilometers', 'radians'}
allUnits = set([param['unitOfMeasure'] for param in params])
for i in np.arange(len(params)):
    if params[i]['unitOfMeasure'] == None:
        continue
    elif params[i]['unitOfMeasure'] == 'kilometers':
        params[i]['unit'] = u.kilometer
    elif params[i]['unitOfMeasure'] == 'km/sec':
        params[i]['unit'] = u.kilometer/u.second
    elif params[i]['unitOfMeasure'] == 'Hz':
        params[i]['unit'] = 1./u.second
    elif params[i]['unitOfMeasure'] == 'degrees':
        params[i]['unit'] = u.degree
    elif params[i]['unitOfMeasure'] == 'rad/min':
        params[i]['unit'] = u.radian/u.minute
    elif params[i]['unitOfMeasure'] == 'Watts':
        params[i]['unit'] = u.joule/u.second
    elif params[i]['unitOfMeasure'] == 'radians':
        params[i]['unit'] = u.radian
    elif params[i]['unitOfMeasure'] == 'units of magnitude':
        params[i]['unit'] = u.mag
    elif params[i]['unitOfMeasure'] == 'decibels':
        params[i]['unit'] = u.decibel
    else:
        print('unhandled unit: ' + str(params[i]['unitOfMeasure']))

