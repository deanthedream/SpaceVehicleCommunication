#grabUDL.py
"""This class defines a function for submitting a URL
and unpacking all data
"""
#call by "import grabUDL.grabUDL as grabUDL" at top of script

import requests, base64, pandas as pd
import astropy.units as u
import numpy as np
from requests.packages.urllib3.exceptions import InsecureRequestWarning

class grabUDL:
    def __init__(self):
        #self.creds = "Basic " + base64.b64encode(b"username:password").decode("ascii")
        #self.creds = "Basic aCharacterStringFromUDLutility=="
        with open("/home/dean/Documents/AFRL2019/myString.txt", 'r') as f:
            myString = f.read().replace("\n","")
        self.creds = "Basic " + myString # myString is in a secured file so this is not pushed into a repo
        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

        self.validQueryParameters = ['obTime','idSensor','maxResults']
        self.validQueryOperations = ['=','=~','=<','=>','in','between']

    def queryUDL(self, url):
        """
        Args:
            url (string) - UDL url query string
        Returns:
            elsetsDataFrame (pandas object) - the pandas data object
            elsetsDataFrame.keys() (list) - list of strings of column names
        """
        result = requests.get(url, headers={'Authorization':self.creds}, verify=False)
        tmp = result.json()
        if hasattr(tmp,'keys'):
            return tmp, tmp.keys()
        else: # it is a list
            elsetsDataFrame = pd.DataFrame(tmp)
            return elsetsDataFrame, elsetsDataFrame.keys()
        #elsetsDataFrame.keys() #Gets all possible columns
        #keepColumns = ['epoch','meanMotion','eccentricity','inclination','meanAnomaly']
        #elsetsArray = elsetsDataFrame[keepColumns].values

    def addAstropyUnits_to_queryHelp(self, queryhelp):
        """Trying to convert these into manipulatable quantities in astropy
        {'units of magnitude', 'Hz', None, 'degrees', 'km/sec', 'rad/min', 'Watts', 'decibels', 'kilometers', 'radians'}
        """
        params = queryhelp['parameters']
        allUnits = set([param['unitOfMeasure'] for param in params])
        for i in np.arange(len(params)):
            if params[i]['unitOfMeasure'] == None:
                continue
            elif params[i]['unitOfMeasure'] == 'kilometers' or \
                 params[i]['unitOfMeasure'] == 'km':
                params[i]['unit'] = u.kilometer
            elif params[i]['unitOfMeasure'] == 'km/sec' or \
                 params[i]['unitOfMeasure'] == 'kilometers per second':
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
            elif params[i]['unitOfMeasure'] == 'degrees per second':
                params[i]['unit'] = u.degree/u.second
            elif params[i]['unitOfMeasure'] == 'seconds':
                params[i]['unit'] = u.second
            elif params[i]['unitOfMeasure'] == 'minutes':
                params[i]['unit'] = u.minute
            elif params[i]['unitOfMeasure'] == 'dBW':
                params[i]['unit'] = u.decibel*u.joule/u.second
            elif params[i]['unitOfMeasure'] == 'dBW/Hz':
                params[i]['unit'] = u.decibel*u.joule
            elif params[i]['unitOfMeasure'] == 'm/s':
                params[i]['unit'] = u.meter/u.second
            elif params[i]['unitOfMeasure'] == 'meters squared':
                params[i]['unit'] = u.meter**2.
            elif params[i]['unitOfMeasure'] == 'revolutions per day':
                #in regards to mean motion so should be julian days
                params[i]['unit'] = 1./u.day
            elif params[i]['unitOfMeasure'] == 'revolutions per day squared':
                #in regards to mean motion so should be julian days
                params[i]['unit'] = 1./u.day**2.
            elif params[i]['unitOfMeasure'] == 'revolutions per day cubed':
                #in regards to mean motion so should be julian days
                params[i]['unit'] = 1./u.day**3.
            elif params[i]['unitOfMeasure'] == 'inverse earth radii':
                params[i]['unit'] = 1./u.earthRad
            else:
                print('unhandled unit: ' + str(params[i]['unitOfMeasure']))
        return params

    def addQueryParameter(self, baseQuery, queryParameter, operation='=', value=None):
        """
        Intended to append a "query parameter" such as idSensor or obTime
        To the UDL query string
        Args:
            baseQuery (string) - the base of the query string to append to
        Returns:
            queryURL (string) - the query URL with the "query parameter" added
        """
        #1 Check if ? or & needs to come next (if no ?, then &)
        if '?' in baseQuery:
            queryURL = '&'
        else:
            queryURL = '?'

        #2 Check if queryParameter is a known parameter and append to URL
        if not queryParameter in self.validQueryParameters:
            print(queryParameter + ' Error: not in validQueryParameters')
        queryURL = queryURL + queryParameter

        #3 operation
        assert operation in self.validQueryOperations, 'Error: operation not in set of operations'
        if operation == 'in':
            # if idSensor in {1,2,3}
            #https://unifieddatalibrary.com/udl/observation?obTime=&idSensor=1%2C2%2C3
            queryURL = queryURL + '='
            for v in value:
                queryURL = queryURL + str(v)
                if not v == value[-1]:
                    queryURL = queryURL + '%2C'
        elif operation == 'between':
            # if idSensor between 1 and 30
            # https://unifieddatalibrary.com/udl/observation?obTime=&idSensor=1..30
            queryURL = queryURL + '=' + str(value[0]) + '..' + str(value[1])
        else:
            queryURL = queryURL + operation + str(value)

        return baseQuery + queryURL