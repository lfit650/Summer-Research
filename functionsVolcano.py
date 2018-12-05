# Function to find temperature at a point y,z underground in a
# Geothermal System.
#import import_Temperature
from ipywidgets import interact, fixed, interactive_output, HBox, Button, VBox, Output, IntSlider, Checkbox, FloatSlider
#import import_Temperature
from ipywidgets import interact, fixed, interactive_output, HBox, Button, VBox, Output, IntSlider, Checkbox, FloatSlider
import scipy as sp
import numpy as np
import numpydoc as npd
from IAPWS.iapws_master.iapws import _utils, IAPWS97
from scipy.interpolate import griddata
import math
import matplotlib.pyplot as plt	
from matplotlib import cm
#import sphinx_rtd_theme as srt

# read temperature data 
x,y,z,T = np.genfromtxt(r'H:\Summer Research\xyz_T.dat', delimiter = ',').T  
	
# get min x slice of data
xmin = np.min(x); tol = 1.e-6
inds = np.where(abs(x-xmin)<tol)
y = y[inds]
z = z[inds]
z = z - np.max(z)
T = T[inds]
	
# interpolate temperatures to finer grid in z direction
ynew = np.unique(y)
	
ymax = 5.e3
dy = 100.
yu = np.unique(y)
i = np.argmin(abs(yu - ymax))		
ynew = list(np.linspace(dy/2.,ymax+dy/2., abs(ymax)/dy+1))+ list(yu[i+1:])
# ye = list(np.linspace(0,ye[i+1], abs(ye[i+1])/dy+1))+ list(ye[i+2:])
	
zmin = -1.4e3
dz = 50.
zu = np.flipud(np.unique(z))
i = np.argmin(abs(zu - zmin))		
znew = list(np.linspace(-dz/2.,zmin-dz/2., abs(zmin)/dz+1))+ list(zu[i+1:])
# ze = list(np.linspace(0,ze[i+1], abs(ze[i+1])/dz+1))+ list(ze[i+2:])
			
[ynew,znew] = np.meshgrid(ynew, znew)
shp = ynew.shape
ynew = ynew.flatten()
znew = znew.flatten()
Tnew = griddata(np.array([y,z]).T,T,np.array([ynew,znew]).T, method='linear', fill_value = np.min(T))
yi,zi,Ti = np.array([ynew,znew,Tnew])

#print(yi)
#print(zi)
#for i in range(61):
#    print(zi[i])
#print(Ti)
#print(len(zi))



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

# print(Phase(2500000, 274))
    
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

    # Chnage xmax for different type of volcanoes (angles).
    xmax = hmax # max radius, distance from centre to edge of volcano. 45 degrees if = hmax
    delh = -(hmax/xmax) * y + hmax
    Pi = Pressure(delh-z)
    return Pi

def testPhase(Pi, Ti):
    if (Phase(Pi,Ti) != "Vapour"):
        return False
    return True
#print(testPhase(2215200, 500))
def get_Gold(hmax):
    # hmax = 
    #x1 = hmax
    # n = 1000
    # rangex = np.linspace(0,x1,num = n, endpoint = False)
    # Pdrop = [0]*n
    # for i in range(len(rangex)):
    #     np.append(Pdrop,pressureDrop(rangex[i],hmax)
    #     #print(Pdrop[i])

    j = 0 
    while (hmax > yi[j]):   # Find values we can use from temperature distribution
        j+=1
    # Initialise arrays    
    rangex = yi[0:j+1]   
    #print(rangex) 
    Pdrop = []
    Pi = []
    Ti = []
    totalVolume = 0
    #count = 0
    for i in range(0, len(rangex)): # Loop through sections in x direction
        Pdrop = np.append(Pdrop,pressureDrop(rangex[i],hmax))
    #print(len(Pdrop))
    for i in range(1, len(rangex)): # Loop through sections in x direction
        count = 0
        for j in range(0,len(yi),61):   #Loop through depths in z direction
            ypoint = rangex[i]
            zpoint = zi[j]
            Pi = np.append(Pi,findPressure(ypoint,zpoint,hmax))
            Ti = np.append(Ti,findTemp(ypoint,zpoint))
            newP = Pi[count] - Pdrop[i]
            #print(Pi)
            #print(newP)
            #print(Ti)
            if (testPhase(newP,Ti[count]+273.15) == False):    # Check if water boils, stop when it does not
                #i = len(rangex)
                volume = math.pi * (-1) *zpoint * ((rangex[i])**2-(rangex[i-1])**2)  # Calculate volume of cylinder
                #print(zpoint)
                #j = len(yi)
                totalVolume += volume
                break
                #print(totalVolume)
            count+=1    
    
    minPorosity = 0.045
    maxPorosity = 0.0966666667
    averagePorosity = 0.07083333333333333
    
    minVolume = totalVolume * minPorosity   # in m^3
    maxVolume = totalVolume * maxPorosity   # in m^3
    averageVolume = totalVolume * averagePorosity   # in m^3 
    
    minGoldConc = 4.78358   # in ppm
    maxGoldConc = 16.33333  # in ppm
    averageGoldConc = 10.55845833333333 # in ppm

    minGold = minVolume * minGoldConc   # in grams
    maxGold = maxVolume * maxGoldConc   # in grams
    averageGold = averageVolume * averageGoldConc   # in grams   

    return averageGold  #in grams


print(get_Gold(1000))
#print(Pdrop) 

# Graphs to display

# Temperature distribution. contour
# Depth where no more boiling occurs for each column. bar/line
# Amount of gold boiled from each column. bar/line
 

 #order temperature. change y and z accordingly?? 
# index = np.argsort(Ti)
# revIndex = index[::-1]
# sortedT = np.sort(Ti)
# revST = sortedT[::-1]
# Sy = copy(yi)
# Sz = copy(zi)
# for i in range(len(Ti)): 
#     Sy[revIndex[i]] = yi[i]
#     Sz[revIndex[i]] = zi[i] 
# #print(sortedT)    
# #print(Sy)
# #print(Sz)
# #print(revST)
# plt.contourf(Sy,Sz,revST)

yi2 = yi.reshape(shp)
plt.contourf(yi.reshape(shp),zi.reshape(shp),Ti.reshape(shp))


# sortedT = np.sort(Ti)
# CS = plt.imshow(sortedT, origin='lower', cmap=cm.jet, extent = [yi.min(),yi.max(),zi.min(),zi.max()], aspect='auto')


# plt.colorbar()
plt.show()