#grabUDL.py
"""This class defines a function for submitting a URL
and unpacking all data
"""
#call by "import grabUDL.grabUDL as grabUDL" at top of script

import requests, base64, pandas as pd


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