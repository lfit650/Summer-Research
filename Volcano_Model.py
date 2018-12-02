# Volcano model
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

def Phase(P,T):
	#Tc = 647 # Kelvin
	#Pc = 22.064 # MPa
	P = P/1000000 # Change to MPa
	#return _utils.getphase(Tc, Pc, T, P, 1, None)
	ph = IAPWS97(P=P,T=T)
	return ph.phase

#print(Phase(2500000, 274))

def findTemp(y,z):
	
	# Function to find temperature at a point y,z underground in a
	# Geothermal System.	
    '''
    Parameter:
        y = (float) distance across from the surface.
        z = (float) distance downwards, below surface.
    Return:
        T = Temperature at point (y,z)
    
    '''

def get_Gold():

'''
Will need, 
Volume beneath volcano.
    Area of volcano (m^2), how deep geothermal waters reach (m). 

Concentration of gold in geothermal waters. (ppm).

Temperature everywhere below area of volcano. 
Pressure evrywhere below area of volcano.
    Can get from temperature distribution??
        This will determine whether water will be boiled or not.
        Volume of water that is boiled. 

Pressure drop will be proportional to height of Volcano above the ground.
    Cone volcanoes, height is always changing        Integrate over dx from max height to 0??

'''

'''
Case 1, Assume cone volcano is symmetrical.
Assume water columns have gerater height when volcano has greater height. (to calculate pressure drop)
Steps of function.

1. Find Pdrop for every column of water for 0 < x < x1 (where width/step size 1 of column is chosen, test with efficiency vs accuracy)
2. Find temperature Ti at points in the column starting at the top and going downwards. (Again step size 2 to be determined from efficiency vs accuracy tests)
3. Calculate Pboil at points using temperature Ti at same point.
4. Find pressure Pi at points in the column starting at the top and going downwards (pressure should be increasing as we get deeper, same step size 2 as temperature)
5. Calculate Pi - Pdrop and check if < Pboil.
6. Once Pi - Pdrop is not < Pboil, then get depth of column to this point andcalculate the volume of water that was boiled into steam from this column. (depth * step size 1 * breadth in z direction??)
7. Loop over all columns until x = x1
8. sum total volume boiled to steam. * 2 because of initial case assumption that volcano was symmetric. 
9. Calculate amount of gold (volume * concentration of gold)

'''

# Return amount of gold boiled in grams = Volume of water boiled * Concentration of Gold	