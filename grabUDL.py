#grabUDL.py
"""This class defines a function for submitting a URL
and unpacking all data
"""
#call by "import grabUDL.grabUDL as grabUDL" at top of script

import requests, base64, pandas as pd
import astropy.units as u
import numpy as np


class grabUDL:
    def __init__(self):
        #self.creds = "Basic " + base64.b64encode(b"username:password").decode("ascii")
        #self.creds = "Basic aCharacterStringFromUDLutility=="
        with open("/home/dean/Documents/AFRL2019/myString.txt", 'r') as f:
            myString = f.read().replace("\n","")
        self.creds = "Basic " + myString # myString is in a secured file so this is not pushed into a repo


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
            elif params[i]['unitOfMeasure'] == 'degrees per second':
                params[i]['unit'] = u.degree/u.second
            elif params[i]['unitOfMeasure'] == 'seconds':
                params[i]['unit'] = u.second
            elif params[i]['unitOfMeasure'] == 'dBW':
                params[i]['unit'] = u.decibel*u.joule/u.second
            elif params[i]['unitOfMeasure'] == 'dBW/Hz':
                params[i]['unit'] = u.decibel*u.joule
            else:
                print('unhandled unit: ' + str(params[i]['unitOfMeasure']))
        return params
