# This file is to model earthquakes causing pressure drops.

import numpy as np
from matplotlib import pyplot as plt
from ipywidgets import interact, fixed, interactive_output, HBox, Button, VBox, Output, IntSlider, Checkbox, FloatSlider
TEXTSIZE = 20
from IPython.display import clear_output
import scipy as sp
import numpydoc as npd
from scipy.interpolate import griddata
from scipy.special import erf as erf
from scipy.optimize import root as root
import math
import matplotlib.pyplot as plt	
from matplotlib import cm
from scipy.integrate import quad as defint

# Initialise all constants

dP0 = -0.1e7    # Boiling threshold Pressure
time = 0.1
a = 1   # Bounds of fracture

#kappa = hydraulic diffusivity
# k = permeability
k = 5.92651154e-14  # m^2
# mu = dynamic viscosity
mu = 0.547e-3   # Ns/m^2
# phi = porosity
phi = 0.070833333333333 #  %
# cf = fluid compressibility
cf = 45.8e-11  # 1/Pa
# cr = rock compressibility
cr = 8.137177997e-10  # 1/Pa
# Or could use total compressibility value ct = 1.125e-9 1/Pa

# dP is pressure drop inside the fracture during the earthquake
dP = 10.e6    
kappa = k/(mu*(phi*(cf+cr)))

def fun(x,time):
    '''
    Returns:
        Function, model + dP0 to be used for root.
    '''
    return (1/2)*dP*(erf((a-x)/(2*np.sqrt(kappa*time))) + erf((a+x)/(2*np.sqrt(kappa*time)))) + dP0

def pressureModel(x,dP,a,kappa,time):
    '''
    Returns:
        Pressure equation.
    '''
    return (1/2)*dP*(erf((a-x)/(2*math.sqrt(kappa*time))) + erf((a+x)/(2*math.sqrt(kappa*time))))

# Code below used to plot hydraulic diffusivity changes effects. 
'''
kappa = k/(mu*(phi*(cf+cr)))
mul=0.1 # multiple of kappa.
# Initialise arrays.
kArray = []
solArray = []
while(mul<=3):  # Loop through diffusivitys upto 3*kappa.
    kappa *= mul
    kArray.append(kappa)    # Append array of different kappa values.
    sol1 = root(fun,x0 = [-1,1], args=(time,))
    solArray.append(sol1.x[1])  # Append solutions/intercepts for different kappa.
    kappa /= mul    # Rest kappa.
    mul += 0.1  # Cycle to next multiple.
# Plotting.    
plt.plot(kArray,solArray,'-r')
#plt.xlabel('Hydraulic Diffusivity',size=TEXTSIZE-2)
#plt.ylabel('Width of Boiling Region', size=TEXTSIZE-2)
#plt.title('Effect of Diffusivity on Boiling Region',size=TEXTSIZE)
#plt.show()
plt.gcf().patch.set_facecolor([66/255.,186/255.,151/255.])
plt.savefig('diffusivityVSregionw.png')#, facecolor = plt.gcf().get_facecolor(), transparent = True)
'''

def plot_EQ(time): 
    '''
    Function showing the curve of pressure drop profiles at different times, where time is a slider in a pyhton notebook.
    Parameters:
        time = time (s)
    Returns:
        Shows plot of pressure model for given time.    
    '''
    # Set up figure window and axes size.
    f = plt.figure(figsize=(12,6))
    ax = plt.axes([0.1,0.1,0.8,0.8])
    # Set plot labels and title.
    ax.set_xlabel('x [m]',size=TEXTSIZE)
    ax.set_ylabel('-P [MPa]', size=TEXTSIZE)
    plt.title('Earthquake Pressure with Varying Times',size=TEXTSIZE)

    # Create x array to input into the model.
    x = np.linspace(-10*a, 10*a, 1000)
    #a = 1.e3
    
    # Below is code to plot for varying times without the slider.
    '''
    pms = []    # Initialise array.
    ts = np.logspace(-1, 3, 10) # Set array of times to cycle through
    for time in ts:
        pm = (1.e-6)*pressureModel(x,-dP,a,kappa,time)
        pms.append(pm) 

    for time,pm in zip(ts,pms):
        plt.plot(x,pm,label='t={:3.2e}s'.format(time))
    '''
    # Comment out if using above code.
    plt.plot(x,pressureModel(x,-dP,a,kappa,time))   # Plot pressure model for given parameters and slider bar of time.
    
    # Plotting threshold line
    plt.axhline((1.e-6)*dP0,linestyle=':',color='k',label='Boiling Threshold Pressure')
    plt.legend(prop={'size': 14})

    # If plotting varying times comment starting here ***

    # Find solution/intercepts between threshold and pressure model using scipy.optimize.root function.
    sol = root(fun,x0 = [-1,1], args=(time,))
    if (sol.success):
        # Create string of solution if intercept was found.
        string = "Solution = ({0:.3f}, {1:.3f})".format(sol.x[0], sol.x[1])
    else:
        string = "No Solution"

    # Plot string as text in figure.
    ax.text(.70, .95, string, position = (0.69,0.05), fontstyle='italic', transform=ax.transAxes, size=TEXTSIZE, color = 'r')
    
    # Set axes limits.
    #ax.set_ylim([-0.8e7,0])
    ax.set_xlim([-10*a,10*a])
    
    # End comment here ***

    #Uncomment if using varying times.
    #plt.gcf().patch.set_facecolor([66/255.,186/255.,151/255.])
    #plt.savefig('EarthquakeVaryingTimes.png', facecolor = plt.gcf().get_facecolor(), transparent = True)
    
    #show plot for notebook with slider.
    plt.show()


def EQslider():
	'''
    Slider function created so time can be a slider if the function is run in a notebook
    '''
	tsldr = IntSlider(value = 16, description='$time$', min=1, max = 35, step=1)
	return VBox([HBox([tsldr]), interactive_output(plot_EQ, {'time':tsldr})])

#plot_EQ(time)    