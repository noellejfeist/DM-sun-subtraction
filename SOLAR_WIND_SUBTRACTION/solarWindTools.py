import numpy as np
from astropy.coordinates import get_sun
from astropy.time import Time
from parReader import Pulsar

#from D. Madison et al. (2018), the model for solar wind
def Z0(timeArr : Time, pulsarObj : Pulsar, returnAngle=False, returnCoords=False):
    #electron density: 1(AU)^2/r^2 where r is the distance from the sun in 3D space
    #integrate over the line of sight from the earth to the pulsar

    #get angular separation
    suns = get_sun(timeArr)
    theta = suns.separation(pulsarObj.SkyCoord()).radian

    #distance from point x (Earth's position in AU) to the sun in the direction of the pulsar
    #1/((x-np.cos(theta))**2 + np.sin(theta)**2)
    #integrate from x=0 to infinity (pulsar is approx. infty away)

    #solution verified with J. Hazboun et al. (2021)
    sintheta = np.sin(theta)
    Z0 = (np.pi - theta)/sintheta * 4.81814e-6 #conversion factor from AU to pc

    #different returns depending if the angular separations and the coordinates should be returned to save computation time
    if (returnAngle and returnCoords):
        return (np.pi - theta)/sintheta * 4.81814e-6, theta, pulsarObj.SkyCoord().gcrs.ra.degree - suns.gcrs.ra.degree, pulsarObj.SkyCoord().gcrs.dec.degree - suns.gcrs.dec.degree
    if (returnAngle):
        return (np.pi - theta)/sintheta * 4.81814e-6, theta
    if (returnCoords):
        return (np.pi - theta)/sintheta * 4.81814e-6, suns.gcrs.ra.degree, suns.gcrs.dec.degree
    return (np.pi - theta)/sintheta * 4.81814e-6

def sinfit(x : float, a : float, b : float): #fit function for Sinusoid
    return a*np.sin(2*np.pi*(x-b)/365)

def expfit(x : float, a : float, b : float, c : float): #fit function for Power Law for the SF
    return a*np.power(x, b-2)+c

def expfit2(x : float, a : float):
    return a*np.power(x, 5/3)

#output sun-subtracted dD values to .dmx file
def to_dmx_file(folder : str, pulsarname : str, p1 : Pulsar, newData : np.ndarray):

    #open file
    output = open(folder + "/" + pulsarname + "_sun-subtracted.dmx", 'w')

    #headings
    output.write("# " + p1.PSR + " dispersion measure variation\n")
    output.write("# Mean DMX value = " + "{:e}".format(np.average(newData)) + "\n")
    output.write("# Uncertainty in average DM = ########\n") #errors in .par and .dmx seemed ambiguous, so this is a place holder
    output.write("# Columns: DMXEP DMX_value DMX_var_err DMXR1 DMXR2")
    if len(p1.DMXF1_ARR) != 0 and len(p1.DMXF2_ARR) != 0: #these are only sometimes included
        output.write(" DMXF1 DMXF2")
    output.write(" DMX_bin\n")

    #print data
    for i in range(len(newData)):
        output.write("{:.4f}".format(p1.DMXEP_ARR[i]) + " " +
            "{: e}".format(newData[i]) + " " +
            "{:.3e}".format(p1.DMX_ERROR_ARR[i]) + " " +
            "{:.4f}".format(p1.DMXR1_ARR[i]) + " " +
            "{:.4f}".format(p1.DMXR2_ARR[i]) + "  ")
        if len(p1.DMXF1_ARR) != 0 and len(p1.DMXF2_ARR) != 0: #only print if the data exists
            output.write("{:.2f}".format(p1.DMXF1_ARR[i]) + " " +
            "{:.2f}".format(p1.DMXF2_ARR[i]) + " ")
        output.write("DX" + str(i + 1).zfill(3) + "\n")

    output.close()