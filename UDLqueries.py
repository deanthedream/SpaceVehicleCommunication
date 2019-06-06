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
dat_help['parameters'] = g.addAstropyUnits_to_queryHelp(dat_help)



#query observation help
url_observationhelp = "https://unifieddatalibrary.com/udl/observation/queryhelp"
dat_observationhelp, col_observationhelp = g.queryUDL(url_observationhelp)
dat_observationhelp['parameters'] = g.addAstropyUnits_to_queryHelp(dat_observationhelp)

#query rfobservation help
url_rfobservationhelp = "https://unifieddatalibrary.com/udl/rfobservation/queryhelp"
dat_rfobservationhelp, col_rfobservationhelp = g.queryUDL(url_rfobservationhelp)
dat_rfobservationhelp['parameters'] = g.addAstropyUnits_to_queryHelp(dat_rfobservationhelp)


