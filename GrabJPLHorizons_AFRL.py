#Grab Data From JPL Horizons
#By: Dean Keithly
#Written: Aug 8,2018

import numpy as np
#from astroquery.jplhorizons import Horizons #vestigal
#import urllib2 as urllib # python2
import urllib.request as urllib
import json
#import re #vestigal
#JPL Horizons query fields defined in ftp://ssd.jpl.nasa.gov/pub/ssd/horizons_batch_example.long

def extractPositions(bodyName, START_TIME, STOP_TIME, STEP_SIZE, OBJ_DATA='YES', MAKE_EPHEM='NO',OUT_UNITS='KM-S',TABLE_TYPE='VECTOR',CENTER='@Sun'):
    """Extracts body positions over the specified time range with the specified step sizes
    Args:
        bodyName (string) - one of the set of body names in bodyNametoNAIFID
        START_TIME (string) - of form 'yyyy-mm-dd' 
        STOP_TIME (string) - of form 'yyyy-mm-dd'
    Returns:
        times (numpy array) - numpy array of times with dims [num times]
        r_body_sun (numpy array) - numpy array of coordinates with dims [num times][3]
    """
    DATA =  {'body_name':bodyName}#Create Data Structure
    DATA['naifID'] = bodyNametoNAIFID(bodyName)#Grab NAIF ID
    html = queryJPLHorizons(DATA['naifID'], START_TIME, STOP_TIME, STEP_SIZE, OBJ_DATA, MAKE_EPHEM, OUT_UNITS, TABLE_TYPE, CENTER)#Query JPL Horizons

    html2 = html.split('\n')
    indexSOE = [i for i in np.arange(len(html2)) if "SOE" in html2[i]][0]
    indexEOE = [i for i in np.arange(len(html2)) if "EOE" in html2[i]][0]
    dataArray = html2[indexSOE+1:indexEOE]
    dataArraySplit = [dataArray[i].split(',') for i in np.arange(len(dataArray))]
    #Data format ['JDTDB','Calendar Date','x','y','z','vx','vy','vz','LT','RG','RR']
    r_body_sun = np.asarray([np.asarray([float(dataArraySplit[i][2]),float(dataArraySplit[i][3]),float(dataArraySplit[i][4])]) for i in np.arange(len(dataArraySplit))]) #Extract r_body_sun
    times = [float(dataArraySplit[i][0]) for i in np.arange(len(dataArraySplit))] #Extract Times

    return times, r_body_sun

def parseEphemeris(html):
    """
    Args:
        html - designed to take in the html returned by a queryJPLHorizons function call
    Return:
        parsedObject
    """
    html2 = html.split('\n')
    return parsedObject

#MIGHT NEED TO CREATE A FUNCTION TO SPLINE THE ABOVE FOR INTERPOLATION PURPOSES...

def extractParams(paramNames,bodyName):
    # """Extracts parameter from JPL Horizons Data Page
    # Args
    #     paramNames (list of strings)
    #     bodyName (string)
    # returns
    #     DATA (struct) contains all the extracted JPL Horizons Data + FlyBy Time Position Vectors
    # """
    DATA =  {'body_name':bodyName}#Create Data Structure
    DATA['naifID'] = bodyNametoNAIFID(bodyName)#Grab NAIF ID
    html = queryJPLHorizons(DATA['naifID'])#Query JPL Horizons

    html = html.replace(",","")#Removing Parameter messing with data extraction (could alternatively change JPL Horizons)
    html = html.replace(" km "," (km) ")#fixing discrepancy in JPL Horizons
    html = html.replace(" kg "," (kg) ")#fixing discrepancy in JPL Horizons

    #Parse the text output for params
    for line in html.split('\n'):#iterate over each line in html
        for param in paramNames:#Iterate over all Parameters
            if isParamInLine(param,line):#Checks if param is in line
                data = extractParamFromLine(param,line)#Extracted data structure for param
                DATA[param] = data

                #### Reduce Parameters with Aphelion, Perihelion, and Mean Values into those parts
                if param in ['Solar Constant (W/m^2)','IR Temperature (K)','Darkside IR Temperature (k)',\
                            'Maximum Planetary IR (W/m^2)','Minimum Planetary IR (W/m^2)']:
                    tmp = DATA[param]#Pull out param to modify
                    if len(tmp) == 3:
                        DATA[param] = {'Aphelion':tmp[0],'Perihelion':tmp[0], 'Mean':tmp[0]}
                if param == 'Solar Constant (W/m^2)' and bodyName == 'MERCURY':
                    print(data)
                    print(DATA)
    return DATA

def bodyNametoNAIFID(bodyName):
    # """Converts from Body Name to naif ID
    # """
    naifIDS = {'MERCURY':199,'VENUS':299,'EARTH':399,'MOON':301,'MARS':499,'PHOBOS':401,'DEIMOS':402,'JUPITER':599,'IO':501,'EUROPA':502,'GANYMEDE':503,\
    'CALLISTO':504,'AMALTHEA':505,'HIMALIA':506,'ELARA':507,'PASIPHAE':508,'SINOPE':509,'LYSITHEA':510,\
    'CARME':511,'ANANKE':512,'LEDA':513,'THEBE':514,'ADRASTEA':515,'METIS':516,'CALLIRRHOE':517,'THEMISTO':518,\
    'MAGACLITE':519,'TAYGETE':520,'CHALDENE':521,'HARPALYKE':522,'KALYKE':523,'IOCASTE':524,'ERINOME':525,\
    'ISONOE':526,'PRAXIDIKE':527,'AUTONOE':528,'THYONE':529,'HERMIPPE':530,'AITNE':531,'EURYDOME':532,\
    'EUANTHE':533,'EUPORIE':534,'ORTHOSIE':535,'SPONDE':536,'KALE':537,'PASITHEE':538,'HEGEMONE':539,\
    'MNEME':540,'AOEDE':541,'THELXINOE':542,'ARCHE':543,'KALLICHORE':544,'HELIKE':545,'CARPO':546,'EUKELADE':547,\
    'CYLLENE':548,'KORE':549,'HERSE':550,'DIA':553,'SATURN':699,'MIMAS':601,'ENCELADUS':602,'TETHYS':603,\
    'DIONE':604,'RHEA':605,'TITAN':606,'HYPERION':607,'IAPETUS':608,'PHOEBE':609,'JANUS':610,'EPIMETHEUS':611,\
    'HELENE':612,'TELESTO':613,'CALYPSO':614,'ATLAS':615,'PROMETHEUS':616,'PANDORA':617,'PAN':618,'YMIR':619,\
    'PAALIAQ':620,'TARVOS':621,'IJIRAQ':622,'SUTTUNGR':623,'KIVIUQ':624,'MUNDILFARI':625,'ALBIORIX':626,\
    'SKATHI':627,'ERRIAPUS':628,'SIARNAQ':629,'THRYMR':630,'NARVI':631,'METHONE':632,'PALLENE':633,\
    'POLYDEUCES':634,'DAPHNIS':635,'AEGIR':636,'BEBHIONN':637,'BERGELMIR':638,'BESTLA':639,'FARBAUTI':640,\
    'FENRIR':641,'FORNJOT':642,'HATI':643,'HYRROKKIN':644,'KARI':645,'LOGE':646,'SKOLL':647,'SURTUR':648,\
    'ANTHE':649,'JARNSAXA':650,'GREIP':651,'TARQEQ':652,'AEGAEON':653,'URANUS':799,'ARIEL':701,'UMBRIEL':702,\
    'TITANIA':703,'OBERON':704,'MIRANDA':705,'CORDELIA':706,'OPHELIA':707,'BIANCA':708,'CRESSIDA':709,\
    'DESDEMONA':710,'JULIET':711,'PORTIA':712,'ROSALIND':713,'BELINDA':714,'PUCK':715,'CALIBAN':716,\
    'SYCORAX':717,'PROSPERO':718,'SETEBOS':719,'STEPHANO':720,'TRINCULO':721,'FRANCISCO':722,'MARGARET':723,\
    'FERDINAND':724,'PERDITA':725,'MAB':726,'CUPID':727,'NEPTUNE':899,'TRITON':801,'NEREID':802,'NAIAD':803,\
    'THALASSA':804,'DESPINA':805,'GALATEA':806,'LARISSA':807,'PROTEUS':808,'HALIMEDE':809,'PSAMATHE':810,\
    'SAO':811,'LAOMEDEIA':812,'NESO':813,'PLUTO':999,'CHARON':901,'NIX':902,'HYDRA':903,'KERBEROS':904,'STYX':905}
    naifID = naifIDS[bodyName]
    return naifID

def queryJPLHorizons(naifID, START_TIME='2000-01-01', STOP_TIME='2000-12-31', STEP_SIZE='15%20d', OBJ_DATA='YES', MAKE_EPHEM='NO',OUT_UNITS='KM-S',TABLE_TYPE='OBSERVER',CENTER='@Sun'):
    # """Queries JPL Horizons for Object Data Page DATA
    # Args:
    #     naifID (int) - naif ID of body
    # Return:
    #     rawHTML (string) - string containing all information responded by url query
    # """
    #List of Input Paramter Names
    inputNames = ['COMMAND','OBJ_DATA','MAKE_EPHEM','TABLE_TYPE','START_TIME','STOP_TIME','STEP_SIZE','QUANTITIES','CSV_FORMAT','OUT_UNITS','CENTER']

    #Struct of Coded Inputs
    inputs = {'COMMAND': str(naifID), 'OBJ_DATA' : OBJ_DATA ,'MAKE_EPHEM': MAKE_EPHEM, 'TABLE_TYPE': TABLE_TYPE, 'START_TIME': START_TIME,
            'STOP_TIME': STOP_TIME, 'STEP_SIZE':STEP_SIZE, 'QUANTITIES': '1,9,20,23,24', 'CSV_FORMAT': 'YES', 'OUT_UNITS':OUT_UNITS, 'CENTER':CENTER}

    #Static Components of URL
    staticString = ["https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND='","'&OBJ_DATA='","'&MAKE_EPHEM='","'&TABLE_TYPE='","'&START_TIME='",\
        "'&STOP_TIME='","'&STEP_SIZE='","'&QUANTITIES='","'&CSV_FORMAT='","'&OUT_UNITS='","'&CENTER='","'"]

    #### Generate URL Command from Inputs and Static Components
    myURL = ''#Initialize Command String
    for i in np.arange(len(inputs.keys())):
        myURL = myURL + staticString[i]#Append Static Command String Component
        myURL = myURL + str(inputs[inputNames[i]])#Append Input
    myURL = myURL + staticString[-1]#Add Last Element
    #print myURL#Display URL for Verification

    #### Submit URL Request
    response = urllib.urlopen(myURL)
    html = response.read().decode('utf-8')

    #print html #Here I print the output from the HTML query
    return html

def isParamInLine(param,line):
    # """Determines if a Parameter Exists in a line
    # Args:
    #     param (string) Parameter name to look for
    #     line (string) string to look in
    # return:
    #     inLine (boolean) True mean it exists in the line
    # """
    #Rule 1: First Word of Param Must Occur Anywhere in line
    #Rule 2: Last Word of Param Must Occur Anywhere in line
    if param.split()[0] in line:
        #Case 1: The entire explicit parameter name is in the line
        if param in line:
            return True
        #Case 2: The last parameter is also in the line
        if param.split()[-1] in line:
            return True
    return False

def stringContainsNumber(string):
    # """Determine if the passed in string is Number
    # Args:
    #     string (string)
    # Return:
    #     isString (boolean) True if string is a number
    # """
    for letter in string:#Check each letter to determine if it is a number
        try:
            float(letter)
            return True
        except Exception:
            continue
    else:
        return False

def stringIsFloat(myFloat):
    # """Checks if string passed in is a float
    # """
    try:
        float(myFloat)
        return True
    except Exception:
        return False

def typeIsFloat(myFloat):
    # """Checks if thing passed in has type float
    # """
    try:
        if type(myFloat) == type(1.0):
            return True
        else:
            return False
    except Exception:
        return False

def extractParamFromLine(param,line): 
    # """Determines if a Parameter Exists in a line
    # Args:
    #     param (string) Parameter name to look for
    #     line (string) string to look in
    # return:
    #     data[{
    #     centerNum
    #     upperNum
    #     loweNum}]
    # """
    #Find Area to start searching For Number From
    paramStart = line.find(param.split()[0])#Start location of first word clump of parameter
    paramEnd = line.find(param.split()[-1]) + len(param.split()[-1])#End location of last word clump of parameter
    line_paramEndToEnd = line[(paramEnd+1):]#string from start of parameter to end of line

    #### Find Multiplier in Param
    multiplier = findMultiplier(line[paramStart:paramEnd])

    #Assume next number glob is a value
    #Next Char must be = OR Number
    data = list()
    if line_paramEndToEnd.split()[0] == '=':#The next Char is =
        #print 'The Next Text Cluster is The Number you want'
        dataString = line_paramEndToEnd.split()[1]
        data.append(parseDataString(dataString,multiplier))
    elif stringContainsNumber(line_paramEndToEnd.split()[0]):#There is nothing between end of param and number
        if any(not stringContainsNumber(word) for word in line_paramEndToEnd.split()): # A word does not contain number
            print('A word in line does not contain number')
            #This means the aphelion, perihelion, mean section of JPL Horizons is seisitive to format and no non-numbers should follow
        else:
            for dataString in line_paramEndToEnd.split():#Iterate over all values 
                data.append(parseDataString(dataString,multiplier))
    return data

def parseDataString(dataString,multiplier=1.):
    # """This method parses a string known to contain data for its +/- values
    # dataString (string) of the follwing forms:
    #     float+-float
    #     float+float-float
    #     float+-float float+-float float+-float
    #     float+float-float float+float-float float+float-float
    #     float float float
    # Args:
    #     dataString
    #     multiplier (float) if the param contains a multiplier, add this to the string
    # """
    if '+-' in dataString:
        splitString = dataString.split('+-')
        assert stringIsFloat(splitString[0])
        assert stringIsFloat(splitString[1])
        centerNum = float(splitString[0])*multiplier
        upperNum = centerNum + float(splitString[1])*multiplier
        lowerNum = centerNum - float(splitString[1])*multiplier
        return [centerNum, upperNum, lowerNum]
    elif '+' in dataString:#Assume'+' and '-' in string
        split1 = dataString.split('+')
        split2 = split1[1].split('-')
        assert stringIsFloat(split1[0])
        assert stringIsFloat(split2[0])
        assert stringIsFloat(split2[1])
        centerNum = float(split1[0])*multiplier
        upperNum = centerNum + float(split2[0])*multiplier
        lowerNum = centerNum - float(split2[1])*multiplier
        return [centerNum, upperNum, lowerNum]
    else:#String only contains a number
        if stringIsFloat(dataString):
            return [float(dataString)*multiplier,None,None]
        else:
            print('string to Handle')
            print(dataString)
            assert False,'there is some parsing error'

def findMultiplier(myActualParamString):
    # """Determines whether there is a x10^YY within observed parameter
    # """
    multiplier = 1.
    for word in myActualParamString.split():#Iterate over all words in string

        if ('x' in word or 'X' in word) and stringContainsNumber(word):#we know a multiplier exists in the parameter
            xInd = word.find('x')#find position of x
            numStr = word[xInd+1:]#get '10^ZZ'
            nums = numStr.split('^')
            multiplier = float(nums[0])**float(nums[1])
    return multiplier

def checkMissingFields(paramNames,bodyDATA,body_names):
    # """Checks Fields for Missing Data
    # Args:
    #     paramNames  (list) - list of strings containing parameter names
    # """
    for lmn in np.arange(len(body_names)):
        body_name = body_names[lmn]
        for param in paramNames:#Iterate over all parameters
            for key in bodyDATA[body_name].keys():#Iterate over all available planetary parameters
                if isParamInLine(param,key):
                    break
            else:
                print('body_name ' + body_name + ' is missing field ' + param)

#Check of body_names is an input to this model
try:
    print(body_names)
    body_names = [str(name) for name in body_names]
except Exception:
    #body_names = ['MERCURY','VENUS','EARTH','MOON','MARS','PHOBOS','DEIMOS','JUPITER']
    body_names = ['MERCURY','VENUS','EARTH','MOON','MARS','PHOBOS','DEIMOS','JUPITER','IO','EUROPA','GANYMEDE',\
        'SATURN','TITAN','NEPTUNE','URANUS','TRITON','PLUTO']
if len([name for name in body_names if not str(name) in [u'','']]) == 0:#IF ANY INPUTS ARE EMPTY STRINGS
    body_names = ['MERCURY','VENUS','EARTH','MOON','MARS','PHOBOS','DEIMOS','JUPITER','IO','EUROPA','GANYMEDE',\
        'SATURN','TITAN','NEPTUNE','URANUS','TRITON','PLUTO']

print(body_names)
#Check of paramNames is an input to this model
try:
    print(paramNames)
    paramNames = [str(param) for param in paramNames]
except Exception:
    paramNames = ['Vol. Mean Radius (km)',\
    'Rot. Rate XX (rad s^-1)',\
    'Solar Constant (W/m^2)',\
    'IR Temperature (K)',\
    'Darkside IR Temperature (k)',\
    'Maximum Planetary IR (W/m^2)',\
    'Minimum Planetary IR (W/m^2)',\
    'Bond Albedo',\
    'Mass XX (kg)',\
    'Surface emissivity',\
    'Geometric Albedo']
if len([param for param in paramNames if not str(param) in [u'','']]) == 0:#IF ANY INPUTS ARE EMPTY STRINGS
    paramNames = ['Vol. Mean Radius (km)',\
    'Rot. Rate XX (rad s^-1)',\
    'Solar Constant (W/m^2)',\
    'IR Temperature (K)',\
    'Darkside IR Temperature (k)',\
    'Maximum Planetary IR (W/m^2)',\
    'Minimum Planetary IR (W/m^2)',\
    'Bond Albedo',\
    'Mass XX (kg)',\
    'Surface emissivity',\
    'Geometric Albedo']
print(paramNames)

bodyDATA= {}#will contain all the information for each body
for lmn in np.arange(len(body_names)):
    body_name = body_names[lmn]#body name we are currently looking at
    bodyDATA[body_name] = extractParams(paramNames,body_name)#Add data for body to structure

    #Uncomment to view missing parameters in JPL Horizons
    #checkMissingFields(paramNames,bodyDATA,[body_name])#Check for parameters missing

#The output of this function is bodyDATA
bodyDATA = json.dumps(bodyDATA)


#GRABBING PLANET POSITIONS AT SPECIFIC TIMES
#synthetic times for testing purposes
# NOTE COULD USE TLIST to get discrete times for query which would need times in MJD
START_TIME = '2019-04-25'
STOP_TIME = '2019-04-30'
STEP_SIZE = '60min'
body_name = 'EARTH'
extractPositions(body_name, START_TIME, STOP_TIME, STEP_SIZE, OBJ_DATA='NO', MAKE_EPHEM='YES', OUT_UNITS='KM-S')