# Function to find temperature at a point y,z underground in a
# Geothermal System.
import import_Temperature
from ipywidgets import interact, fixed, interactive_output, HBox, Button, VBox, Output, IntSlider, Checkbox, FloatSlider
import import_Temperature
from ipywidgets import interact, fixed, interactive_output, HBox, Button, VBox, Output, IntSlider, Checkbox, FloatSlider
import scipy as sp
import numpy as np
import numpydoc as npd
from IAPWS.iapws_master.iapws import _utils, IAPWS97
#import sphinx_rtd_theme as srt

def Pressure(z): 
	'''
	P: Pressure beneath volcano.
	rho: Density of water.
	g: Gravity.
	z: Height of water beneath the volcano.

	'''
	rho = 1000
	g = 9.81
	P = rho*g*z
	return P

def pressureDrop(y,hmax):
	# hmax = 				# max height of volcano
    xmax = hmax				# max radius, distance from centre to edge of volcano. 
    delh = -(hmax/xmax) * y + hmax
    Pdrop = Pressure(delh)
    return Pdrop  

def Phase(P,T):
	#Tc = 647 # Kelvin
	#Pc = 22.064 # MPa
	P = P/1000000 # Change to MPa
	#return _utils.getphase(Tc, Pc, T, P, 1, None)
	ph = IAPWS97(P=P,T=T)
	return ph.phase

#print(Phase(2500000, 274))
    
def findTemp(y,z):
    '''
    Parameter:
        y = (float) distance across from the surface.
        z = (float) distance downwards, below surface.
    Return:
        T = Temperature at point (y,z)
    
    '''
    i = 0
    while (zi[i] != z):
        i+=1
    j = i
    while (yi[j] != y and j<i+61):
        j+=1
    T = Ti[j]
    return T

def findPressure(y,z,hmax):
    # hmax = max height of volcano
    # xmax =  max radius, distance from centre to edge of volcano. 
    delh = -(hmax/xmax) * y + hmax
    Pi = Pressure(delh-z)
    return Pi

def testPhase(Pi, Ti):
    if (Phase(Pi,Ti) != "vapour"):
        return False
    return True

def get_Gold(hmax):
    # hmax = 
    #x1 = hmax
    '''
    n = 1000
    rangex = np.linspace(0,x1,num = n, endpoint = False)
    Pdrop = [0]*n
    for i in range(len(rangex)):
        np.append(Pdrop,pressureDrop(rangex[i],hmax)
        #print(Pdrop[i])
    '''
    j = 0 
    while (hmax > yi[j]):   # Find values we can use from temperature distribution
        j+=1
    # Initialise arrays    
    rangex = yi[0:j+1]    
    Pdrop = [0]*len(rangex)
    Pi = [0]*len(rangex)
    Ti = [0]*len(rangex)
    totalVolume = 0
    #count = 0
    for i in range(1, len(rangex)): # Loop through sections in x direction
        count = 0
        np.append(Pdrop,pressureDrop(rangex[i],hmax))
        for j in range(0,len(yi),61):   #Loop through depths in z direction
            ypoint = rangex[i]
            zpoint = zi[j]
            np.append(Pi,findPressure(ypoint,zpoint,hmax))
            np.append(Ti,findTemp(ypoint,zpoint))
            newP = Pi[count] - Pdrop[count]
            if (testPhase(newP,Ti[count]) == False):    #Check if water boils, stop when it does not
                #i = len(rangex)
                volume = math.Pi * (rangex[i]-rangex[i-1])**2 * zpoint  # Calculate volume of cylinder
                j = len(yi)
                totalVolume += volume
            count+=1    


#get_Gold(2797)
#print(Pdrop)