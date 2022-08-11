from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

#finds data with appropriate heading, adds to list called varName
def extractData(line, heading, varName, argument, exactHeading=False, changeToFloat=True):
    #normally only looks for "heading" in the variable heading, but if exactHeading, they must match exactly

    lineSplit = line.split()

    if (exactHeading):
        if lineSplit[0] != heading:
            return
        data = lineSplit[argument]
        if (changeToFloat): 
            data = float(data.replace('D','E')) #reformat to change to a float
        varName.append(data)
        return
    
    try:
        lineSplit[0].index(heading)
        data = lineSplit[argument]
        if (changeToFloat): 
            data = float(data.replace('D','E')) #reformat to change to a float
        varName.append(data)
        return
    except:
        return

def singleVarExtractData(line, heading, varName, argument, changeToFloat=True): #similar as before, modified for vars instead of lists
    lineSplit = line.split()
    
    try:
        a = lineSplit[0].index(heading)
        if a != 0:
            return varName
        data = lineSplit[argument]
        if (changeToFloat): 
            data = float(data.replace('D','E')) #reformat to change to a float
        return data
    except:
        return varName

#class to store .par information
class Pulsar:

    def __init__(self, filename): #constructor

        file = open(filename, 'r')

        self.Init() #initialize all variables to None and all lists to []

        for line in file: #iterate through all the DM measurements, and store appropriate variables

            self.PSR = singleVarExtractData(line, "PSR", self.PSR, 1, changeToFloat=False)
            self.LAMBDA = singleVarExtractData(line, "LAMBDA", self.LAMBDA, 1)
            self.LAMBDA_ERROR = singleVarExtractData(line, "LAMBDA", self.LAMBDA_ERROR, 3)
            self.BETA = singleVarExtractData(line, "BETA", self.BETA, 1)
            self.BETA_ERROR = singleVarExtractData(line, "BETA", self.BETA_ERROR, 3)
            self.ELONG = singleVarExtractData(line, "ELONG", self.ELONG, 1)
            self.ELONG_ERROR = singleVarExtractData(line, "ELONG", self.ELONG_ERROR, 3)
            self.ELAT = singleVarExtractData(line, "ELAT", self.ELAT, 1)
            self.ELAT_ERROR = singleVarExtractData(line, "ELAT", self.ELAT_ERROR, 3)
            self.PMLAMBDA = singleVarExtractData(line, "PMLAMBDA", self.PMLAMBDA, 1)
            self.PMLAMBDA_ERROR = singleVarExtractData(line, "PMLAMBDA", self.PMLAMBDA_ERROR, 3)
            self.PMBETA = singleVarExtractData(line, "PMBETA", self.PMBETA, 1)
            self.PMBETA_ERROR = singleVarExtractData(line, "PMBETA", self.PMBETA_ERROR, 3)
            self.PX = singleVarExtractData(line, "PX", self.PX , 1)
            self.PX_ERROR = singleVarExtractData(line, "PX", self.PX_ERROR, 3)
            self.POSEPOCH = singleVarExtractData(line, "POSEPOCH", self.POSEPOCH, 1)
            self.F0 = singleVarExtractData(line, "F0", self.F0, 1)
            self.F0_ERROR = singleVarExtractData(line, "F0", self.F0_ERROR, 3)
            self.F1 = singleVarExtractData(line, "F1", self.F1, 1)
            self.F1_ERROR = singleVarExtractData(line, "F1", self.F1_ERROR, 3)
            self.PEPOCH = singleVarExtractData(line, "PEPOCH", self.PEPOCH, 1)
            self.START = singleVarExtractData(line, "START", self.START, 1)
            self.FINISH = singleVarExtractData(line, "FINSIH", self.FINSIH, 1)
            self.DM = singleVarExtractData(line, "DM", self.DM, 1)
            self.DMX = singleVarExtractData(line, "DMX", self.DMX, 1)
            extractData(line, "DMX_", self.DMX_ARR, 1)
            extractData(line, "DMX_", self.DMX_ERROR_ARR, 3)
            extractData(line, "DMXEP_", self.DMXEP_ARR, 1)
            extractData(line, "DMXR1_", self.DMXR1_ARR, 1)
            extractData(line, "DMXR2_", self.DMXR2_ARR, 1)
            extractData(line, "DMXF1_", self.DMXF1_ARR, 1)
            extractData(line, "DMXF2_", self.DMXF2_ARR, 1)
            self.FD1 = singleVarExtractData(line, "FD1", self.FD1, 1)
            self.FD1_ERROR = singleVarExtractData(line, "FD1", self.FD1_ERROR, 3)
            self.FD2 = singleVarExtractData(line, "FD2", self.FD2, 1)
            self.FD2_ERROR = singleVarExtractData(line, "FD2", self.FD2_ERROR, 3)
            self.FD3 = singleVarExtractData(line, "FD3", self.FD3, 1)
            self.FD3_ERROR = singleVarExtractData(line, "FD3", self.FD3_ERROR, 3)
            self.SOLARN0 = singleVarExtractData(line, "SOLARN0", self.SOLARN0, 1)
            self.EPHEM = singleVarExtractData(line, "EPHEM", self.EPHEM, 1, changeToFloat=False)
            self.ECL = singleVarExtractData(line, "ECL", self.ECL, 1, changeToFloat=False)
            self.CLK = singleVarExtractData(line, "CLK", self.CLK, 1, changeToFloat=False)
            self.UNITS = singleVarExtractData(line, "UNITS", self.UNITS, 1, changeToFloat=False)
            self.TIMEEPH = singleVarExtractData(line, "TIMEEPH", self.TIMEEPH, 1, changeToFloat=False)
            self.T2CMETHOD = singleVarExtractData(line, "T2CMETHOD", self.T2CMETHOD, 1, changeToFloat=False)
            self.EPHEM = singleVarExtractData(line, "EPHEM", self.EPHEM, 1, changeToFloat=False)
            self.CORRECT_TROPOSPHERE = singleVarExtractData(line, "CORRECT_TROPOSPHERE", self.CORRECT_TROPOSPHERE, 1, changeToFloat=False)
            self.PLANET_SHAPIRO = singleVarExtractData(line, "PLANET_SHAPIRO", self.PLANET_SHAPIRO, 1, changeToFloat=False)
            self.DILATEFREQ = singleVarExtractData(line, "DILATEFREQ", self.DILATEFREQ, 1, changeToFloat=False)
            self.NTOA = singleVarExtractData(line, "NTOA", self.NTOA, 1, changeToFloat=False)
            self.TRES = singleVarExtractData(line, "TRES", self.TRES, 1)
            self.TZRMJD = singleVarExtractData(line, "TZRMJD", self.TZRMJD, 1)
            self.TZRFRQ = singleVarExtractData(line, "TZRFRQ", self.TZRFRQ, 1)
            self.TZRSITE = singleVarExtractData(line, "TZRSITE", self.TZRSITE, 1, changeToFloat=False)
            self.MODE = singleVarExtractData(line, "MODE", self.MODE, 1, changeToFloat=False)
            self.NITS = singleVarExtractData(line, "NITS", self.NITS, 1, changeToFloat=False)
            self.BINARY = singleVarExtractData(line, "BINARY", self.BINARY, 1, changeToFloat=False)
            self.A1 = singleVarExtractData(line, "A1", self.A1, 1)
            self.A2_ERROR = singleVarExtractData(line, "A1", self.A1_ERROR, 3)
            self.PB = singleVarExtractData(line, "PB", self.PB, 1)
            self.PB_ERROR = singleVarExtractData(line, "PB", self.PB_ERROR, 3)
            self.TASC = singleVarExtractData(line, "TASC", self.TASC, 1)
            self.TASC_ERROR = singleVarExtractData(line, "TASC", self.TASC_ERROR, 3)
            self.EPS1 = singleVarExtractData(line, "EPS1", self.EPS1, 1)
            self.EPS1_ERROR = singleVarExtractData(line, "EPS1", self.EPS1_ERROR, 3)
            self.SINI = singleVarExtractData(line, "SINI", self.SINI, 1)
            self.SINI_ERROR = singleVarExtractData(line, "SINI", self.SINI_ERROR, 3)
            self.M2 = singleVarExtractData(line, "M2", self.M2, 1)
            self.M2_ERROR = singleVarExtractData(line, "M2", self.M2_ERROR, 3)
            self.RNAMP = singleVarExtractData(line, "RNAMP", self.RNAMP, 1)
            self.RNIDX = singleVarExtractData(line, "RNIDX", self.RNIDX, 1)

        file.close()

        #convert to numpy arrays
        self.DMX_ARR = np.asarray(self.DMX_ARR)
        self.DMX_ERROR_ARR = np.asarray(self.DMX_ERROR_ARR)
        self.DMXEP_ARR = np.asarray(self.DMXEP_ARR)
        self.DMXR1_ARR = np.asarray(self.DMXR1_ARR)
        self.DMXR2_ARR = np.asarray(self.DMXR2_ARR)
        self.DMXF1_ARR = np.asarray(self.DMXF1_ARR)
        self.DMXF2_ARR = np.asarray(self.DMXF2_ARR)

        #two possible names for coordinates, set them to match
        if self.ELONG == None:
            self.ELONG = self.LAMBDA
            self.ELAT = self.BETA
        elif self.LAMBDA == None:
            self.LAMBDA = self.ELONG
            self.BETA = self.ELAT

        #create DMXEP_ARR if it does not exist, as this is what I use in solarWindPlots.py
        if len(self.DMXEP_ARR) == 0:
            self.DMXEP_ARR = (self.DMXR1_ARR + self.DMXR2_ARR)/2
        

    def SkyCoord(self): #returns SkyCoord object for pulsar's location
        return SkyCoord(lon=self.LAMBDA*u.degree, lat=self.BETA*u.degree, frame="barycentricmeanecliptic")

    def Init(self): #set all vars to None and all lists to []
        self.PSR = None
        self.LAMBDA = None
        self.LAMBDA_ERROR = None
        self.ELONG = None
        self.ELONG_ERROR = None
        self.BETA = None
        self.BETA_ERROR = None
        self.ELAT = None
        self.ELAT_ERROR = None
        self.PMLAMBDA = None
        self.PMLAMBDA_ERROR = None
        self.PMBETA = None
        self.PMBETA_ERROR = None
        self.PX = None
        self.PX_ERROR = None
        self.POSEPOCH = None
        self.F0 = None
        self.F0_ERROR = None
        self.F1 = None
        self.F1_ERROR = None
        self.PEPOCH = None
        self.START = None
        self.FINSIH = None
        self.DM = None
        self.DMX = None
        self.FD1 = None
        self.FD1_ERROR = None
        self.FD2 = None
        self.FD2_ERROR = None
        self.FD3 = None
        self.FD3_ERROR = None
        self.SOLARN0 = None
        self.EPHEM = None
        self.ECL = None
        self.CLK = None
        self.UNITS = None
        self.TIMEEPH = None
        self.T2CMETHOD = None
        self.EPHEM = None
        self.CORRECT_TROPOSPHERE = None
        self.PLANET_SHAPIRO = None
        self.DILATEFREQ = None
        self.NTOA = None
        self.TRES = None
        self.TZRMJD = None
        self.TZRFRQ = None
        self.TZRSITE = None
        self.MODE = None
        self.NITS = None
        self.BINARY = None
        self.A1 = None
        self.A1_ERROR = None
        self.PB = None
        self.PB_ERROR = None
        self.TASC = None
        self.TASC_ERROR = None
        self.EPS1 = None
        self.EPS1_ERROR = None
        self.SINI = None
        self.SINI_ERROR = None
        self.M2 = None
        self.M2_ERROR = None
        self.RNAMP = None
        self.RNIDX = None
        self.DMX_ARR = []
        self.DMX_ERROR_ARR = []
        self.DMXEP_ARR = []
        self.DMXR1_ARR = []
        self.DMXR2_ARR = []
        self.DMXF1_ARR = []
        self.DMXF2_ARR = []