"""
The purpose of this script is to take a spacecraft TLE
and correlate a specific ground observation's photometric and astrometric
uncertainties to quantify the spacecraft size
"""

import numpy as np
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv

# For a given TLE, try running 
# sgp4 using from sgp4.io import twoline2rv
#and pyorbital using from pyorbital import tlefile

with open("/home/dean/Documents/AFRL2019/tle-new.txt", 'r') as f:
    # read all TLE from file and remove final line caused by split
    lines = f.read().split('\n')[:-1] 

assert np.mod(len(lines),3) == 0, 'the number of lines for each TLE is not 3'
numTLE = len(lines)/3 # there should be 1 TLE for every 3 lines

satData = {} # stores all parsed satellite data
satNums = list() # stores all satellite identification numbers
for i in np.arange(numTLE, dtype='int'):
    line = lines[i*3]
    satName = line.replace(" ","")
    assert len(line) == 24, 'name line must be 24 chars long'
    TLE1 = lines[i*3+1] # first line of TLE
    TLE2 = lines[i*3+2] # second line of TLE

    satOBJ = twoline2rv(TLE1, TLE2, wgs72)
    satNums.append(satOBJ.satnum)
    satData[satOBJ.satnum] = {
    "satName":satName,
    "satelliteOBJ":satOBJ,
    "TLE0":line,
    "TLE1":TLE1,
    "TLE2":TLE2
    }

