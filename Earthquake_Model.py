import numpy as np
from matplotlib import pyplot as plt
from ipywidgets import interact, fixed, interactive_output, HBox, Button, VBox, Output, IntSlider, Checkbox, FloatSlider
TEXTSIZE = 16
from IPython.display import clear_output
import time
import scipy as sp
import numpydoc as npd
from scipy.interpolate import griddata
import math
import matplotlib.pyplot as plt	
from matplotlib import cm

x,y,z,T = np.genfromtxt(r'H:\Summer Research\xyz_T.dat', delimiter = ',').T  
	
# get min x slice of data
xmin = np.min(x); tol = 1.e-6
inds = np.where(abs(x-xmin)<tol)
y = y[inds]
z = z[inds]
z = z - np.max(z)
P = P[inds]
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
Pnew = griddata(np.array([y,z]).T,T,np.array([ynew,znew]).T, method='linear', fill_value = np.min(T))
yi,zi,Pi = np.array([ynew,znew,Tnew])

def findPressure(y,z):

    ### To find Pressure dP at a  ###

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
    P = Pi[j]
    return dP


def plot_EQ(time): # a boundary as a parameter???

    f = plt.figure(figsize=(12,6))
    ax = plt.axes([0.1,0.1,0.8,0.8])
    ax.set_xlabel('x []',size=TEXTSIZE)
    ax.set_ylabel('P []', size=TEXTSIZE)

    #kappa = hydraulic diffusivity
    # k = permeability
    k = 
    # mu = dynamic viscosity
    mu = 
    # phi = porosity
    phi = 
    # cf = fluid compressibility
    cf = 
    # cr = rock compressibility
    cr = 
    # dP instead of T0 = pressure drop inside the fracture during the earthquake
    dP = 
    
    a = #bounds
    kappa = k/(mu*(phi*cf+cr)
    
    x = np.linspace(-a, a, 1000)
    ax.plot(x, ((1/2)*dP*(math.erf((a-x)/(2*math.sqrt(kappa*time)))) + (math.erf((a+x)/(2*math.sqrt(kappa*time))))))

    plt.show()
    ''' 
    At temperature T, find Pboiling (Research or thermodynamic properties library IAPWS)
    Find if Pdrop satisfies equation: P0-Pdrop < Pboiling.
    From this we find volume of water boiled.
    Therefore amount of gold produced.

    '''

def EQslider():
	
	tsldr = IntSlider(value = 50, description='$time$', min=10, max = 100, step=10)
	return VBox([HBox([tsldr]), interactive_output(plot_EQ, {'time':tsldr})])

plot_EQ(20)    