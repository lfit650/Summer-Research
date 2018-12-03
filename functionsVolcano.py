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

def pressureDrop(y):
	# hmax = 				# max height of volcano
    # xmax = 				# max radius, distance from centre to edge of volcano. 
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

def findPressure(y,z):
    # hmax = max height of volcano
    # xmax =  max radius, distance from centre to edge of volcano. 
    delh = -(hmax/xmax) * y + hmax
    Pi = Pressure(delh-z)
    return Pi

def testPhase(Pi, Ti)
    if (Phase(Pi,Ti) == "vapour"):
        return True
    return False

def get_Gold():


  



    






def tempDistribution():
    '''
    Display Temperature distribution graph.
    Sliders for points.
    
    '''
    ysldr = IntSlider(value = 50, description='$yi$', min=50, max = 9750, step=50)
    zsldr = IntSlider(value = -3800, description='$zi$', min=-3800, max = -25, step=50)

#tempDistribution()    

