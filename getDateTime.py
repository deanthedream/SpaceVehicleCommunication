# get Date Time

from datetime import datetime, timedelta


def getULD_UTC(timeDel=0.):
    """
    Args:
        timeDel (float) - the time delta in days. If you want 5 days ago
            timeDel=-5.
    Returns:
        UDL_UTCtime (string) - current UTC time in UDL format
    """
    UTC = datetime.utcnow() + timedelta(days=timeDel)
    UDL_UTCtime = UTC.strftime("%Y-%m-%dT%H:%M:%S.0000Z") # because UDL has a random T, decimals, and Z
    #print(UDL_UTCtime)
    return UDL_UTCtime


#Create a function to parse UDL times from observations


def formUDLtimeString(UTCtime):
    """ Given the UTC time input below, get the obTime in the URL
    # MJD input to utility 1240.65165545
    # UTC time copied from utility 4710-04-18T03:38:23.030879Z
    #https://unifieddatalibrary.com/udl/observation?id=&obTime=>4710-04-18T03%3A38%3A23.030879Z
    Args:
        UTCtime (string) - a string of UTC time the the format in this header
    Returns:
        UTCout (string) - the time string in the UDL format in the header
    """
    goal = UTCtime.replace(":","%3A")
    return goal

