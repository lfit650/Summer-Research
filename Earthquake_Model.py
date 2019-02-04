import numpy as np
from matplotlib import pyplot as plt
from ipywidgets import interact, fixed, interactive_output, HBox, Button, VBox, Output, IntSlider, Checkbox, FloatSlider
TEXTSIZE = 16
from IPython.display import clear_output
#import time
import scipy as sp
import numpydoc as npd
from scipy.interpolate import griddata
from scipy.special import erf as erf
from scipy.optimize import root as root
import math
import matplotlib.pyplot as plt	
from matplotlib import cm

x,y,z,P = np.genfromtxt(r'H:\Summer Research\31kBlockEWA_BH_final_xyz_P.dat', delimiter = ',').T 
	
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
Pnew = griddata(np.array([y,z]).T,P,np.array([ynew,znew]).T, method='linear', fill_value = np.min(P))
yi,zi,Pi = np.array([ynew,znew,Pnew])


#print(yi)

def findPressure(y,z):

    ### To find Pressure dP at a  ###

    '''
    Parameter:
        y = (float) distance across from the surface.
        z = (float) distance downwards, below surface.
    Return:
        P = Pressure at point (y,z)
    
    '''
    i = 0
    while (zi[i] != z):
        i+=1
    j = i
    while (yi[j] != y and j<i+61):
        j+=1
    P = Pi[j]

    
    return P    # As pressure drops to zero when fault opens
dP0 = -0.1e7
#time = 0.1
a = 1 
#a = math.inf #bounds/infinite

#kappa = hydraulic diffusivity
# k = permeability
k = 5.92651154e-14  # m^2
# mu = dynamic viscosity: at 50 deg Clesius what's temperature?
mu = 0.547e-3   # Ns/m^2
# phi = porosity
phi = 0.070833333333333 #  %
# cf = fluid compressibility: what's temperature?
cf = 45.8e-11  # 1/Pa
# cr = rock compressibility
cr = 8.137177997e-10  # 1/Pa
# Or could use total compressibility value ct = 1.125e-9 1/Pa

# dP instead of T0 = pressure drop inside the fracture during the earthquake
# Cycle through full depth to get all pressures??
# Get largest pressure drop.
'''
dP = Pi[0]
#dP = 10.
for i in range (1,len(Pi)):
    if (Pi[i]>dP):
        dP = Pi[i]
'''
dP = 10.e6
#dP = 101325 # Atmopspheric Pressure    
kappa = k/(mu*(phi*(cf+cr)))
def fun(x,time):
    return (1/2)*dP*(erf((a-x)/(2*np.sqrt(kappa*time))) + erf((a+x)/(2*np.sqrt(kappa*time)))) + dP0

def pressureModel(x,dP,a,kappa,time):
    return (1/2)*dP*(erf((a-x)/(2*math.sqrt(kappa*time))) + erf((a+x)/(2*math.sqrt(kappa*time))))

def plot_EQ(time): # a boundary as a parameter???

    f = plt.figure(figsize=(12,6))
    ax = plt.axes([0.1,0.1,0.8,0.8])
    ax.set_xlabel('x []',size=TEXTSIZE)
    ax.set_ylabel('P []', size=TEXTSIZE)

    #print(kappa)
    #fig = plt.figure()
    #ax = plt.axes()

    x = np.linspace(-10*a, 10*a, 1000)
    #a = 1.e3
    '''
    pms = []
    ts = np.logspace(-1, 3, 10)
    for time in ts:
        pm = pressureModel(x,-dP,a,kappa,time)
        pms.append(pm)
        # CHECK erf cannot take array x as input, only scalar number. sp.special works, graph looks incorrect. 
    
    #for i in range (0,len(x)):
    for time,pm in zip(ts,pms):
        plt.plot(x,pm,label='t={:3.2e}'.format(time))
    '''
    plt.plot(x,pressureModel(x,-dP,a,kappa,time))
    plt.axhline(dP0,linestyle=':',color='k')
    #plt.legend()
    sol = root(fun,x0 = [-1,1], args=(time,))
    if (sol.success):
        string = "Solution = ({0:.3f}, {1:.3f})".format(sol.x[0], sol.x[1])
    else:
        string = "No Solution"

    ax.text(.70, .95, string, position = (0.69,0.05), fontstyle='italic', transform=ax.transAxes, size=TEXTSIZE, color = 'r')
    ax.set_ylim([-0.8e7,0])
    #print(sol.x)
    #print(sol.success)
    #plt.savefig('EarthquakeVaryingTime.png')
    plt.show()


def EQslider():
	
	tsldr = IntSlider(value = 16, description='$time$', min=1, max = 35, step=1)
	return VBox([HBox([tsldr]), interactive_output(plot_EQ, {'time':tsldr})])



#plot_EQ(time)    