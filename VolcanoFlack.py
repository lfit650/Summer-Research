from ipywidgets import interact, fixed, interactive_output, HBox, Button, VBox, Output, IntSlider, Checkbox, FloatSlider
import scipy as sp
import numpy as np
import numpydoc as npd
from IAPWS.iapws_master.iapws import _utils, IAPWS97, iapws97
from scipy.interpolate import griddata
import math
import matplotlib.pyplot as plt	
from matplotlib import cm
from scipy.optimize import curve_fit
from matplotlib import pylab
#import sphinx_rtd_theme as srt

'''
Import Pressure and Temperature data
Returns:
    4 1d arrays
    yi = array of distances in the horizontal direction (m)
    zi = distances of depth (m)
    Ti = array of temperatures at points corresponding to yi and zi (degrees Celsius)
    Pi = array of pressuress at points corresponding to yi and zi (Pa)
'''    

# read and import temperature and pressure data 
x,y,z,P = np.genfromtxt(r'H:\Summer Research\31kBlockEWA_BH_final_xyz_P.dat', delimiter = ',').T  
x,y,z,T = np.genfromtxt(r'H:\Summer Research\xyz_T.dat', delimiter = ',').T  
	
# get min x slice of data
xmin = np.min(x); tol = 1.e-6
inds = np.where(abs(x-xmin)<tol)
y = y[inds]
z = z[inds]
z = z - np.max(z)
T = T[inds]
P = P[inds]	
# interpolate temperatures/pressures to finer grid in z direction
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
Pnew = griddata(np.array([y,z]).T,P,np.array([ynew,znew]).T, method='linear', fill_value = np.min(P))
yi,zi,Pi = np.array([ynew,znew,Pnew])

def Pressure(z): 
	'''
	rho = Density of water.
	g = Gravity.
    Parameters:
    	z = Height of water column (m)
    returns:
        P = Pressure (Pa)
	'''
	rho = 1000
	g = 9.81
	P = rho*g*z
	return P

def pressureDrop(y,xmax,hmax):
    '''
    Parameters:
        y = Horizontal distance away from centre (m)
        xmax = Max length of volcano from centre to edge (m)
        hmax = Max height of volcano (m)
    Returns:
        Pressure drop from volcano collapse at this y point (Pa)    
    '''
    if (y>=(3*xmax/5)):
        delh = -(hmax/xmax)*y + hmax
    # else statement needed for flank collapse as not the whole volcano collapses
    # so we calculate the difference in height of original volcano and collapsed volcano.    
    else:
        delh = -(hmax/xmax)*y + hmax - (-(hmax/xmax)*y + 3*hmax/5)
    # Calculate Pressure Drop for new calculated height of collapse.
    Pdrop = Pressure(delh)

    return Pdrop

    
def findTemp(y,z):
    '''
    Parameter:
        y = (float) horizontal distance away from centre (m)
        z = (float) depth below surface (m)
    Return:
        T = Temperature at point (y,z) (degrees Celsius)
    '''
    # Find height in array of heights in distribution.
    i = 0
    while (zi[i] != z):
        i+=1
    j = i
    # Find y point in array of distances in distribution.
    while (yi[j] != y and j<i+61):
        j+=1
    
    # Get Temperature corresponding to this point.
    T = Ti[j]

    return T

'''
Solubility data due to change in temperature only at different pH levels given
as 2 arrays:
    xTemp = Array of temperatures (degrees Celsius)
    yConc = Array of concentrations of gold (ppm)
'''
# At 20 degrees C pH of 7.45.
xTemp = np.array([165, 200, 213.333333, 230, 242, 252, 268.5, 279.8, 290, 296])
yConc = np.array([5, 10, 15, 20, 25, 30, 40, 50, 60, 67])
    
# At 20 degrees C pH of 9.55.
# xTemp = np.array([164.9, 190, 240, 260, 268, 280, 284.5, 290, 298.3, 300])
# yConc = np.array([2.5, 3, 5, 8.5, 10, 13.5, 15, 17, 20, 21.25])

def exp_Fit(x, a, b, c):
    '''
    Returns:
        exponential function.
    '''
    
    return a*np.exp(b*x)+c

# curve_fit function used to fit given xTemp and yConc arrays to an exponential function,
# where popt is an array of optimal parameters.
popt, pcov = curve_fit(exp_Fit, xTemp, yConc, p0=(1, 1e-6, 1))

# Code below used to plot fitted exponential curve 
'''
xxTemp = np.linspace(0, 320, 2501)
yyConc = exp_Fit(xxTemp, *popt)

plt.plot(xTemp,yConc,'o', xxTemp, yyConc)
plt.title('Gold Solubility Changes due to Temperature')
plt.xlabel('Temp. (\xb0C)')
plt.ylabel('Au Conc. (ppm)')
plt.savefig('expFit.png')
'''

def get_Gold(xmax,hmax):
    '''
    Function to give the amount of gold that is deposited from the volcano flank collapse.
    Parameters:
        xmax = Max length of volcano from centre to edge (m)
        hmax = Max height of volcano (m)
    Returns:
        averageGold = Average amount of gold deposited (tonnes)
    '''
    # Researched values for Porosity
    minPorosity = 0.045
    maxPorosity = 0.0966666667
    averagePorosity = 0.07083333333333333
    
    # Initial values of gold.
    minGold = 0
    maxGold = 0
    averageGold = 0

    j = 0
    while (xmax >= yi[j]):   # Cycle through horizontal distances to find values we can use from
        j+=1                 # temperature distribution, that is values that are beneath the volcano.

    # All yi values now in array rangex   
    rangex = yi[0:j]
    
    # Calculate Pressure drop at all distances for given volcano dimensions 
    Pdrop = [0]*2501    # Initialise to correct length.

    for i in range(0, len(rangex)): # Loop through sections in x direction (yi)
            j=i
            Pdrop[j] = (1.e-6)*pressureDrop(rangex[i],xmax,hmax)
            # Calculate pressure drop for all depths in distribution at this x point.
            # depths are 61 indices away in array for same x points.  
            while (j+61<len(yi)):
                j+=61
                Pdrop[j] = (1.e-6)*pressureDrop(rangex[i],xmax,hmax)
    
    # Make array as easier to use later on.
    PdropArray = np.asarray(Pdrop)
    
    # Calculate pressure required to boil at all temperatures through out distribution.
    # This pressure is saturation pressure.
    Parray = [] # Initialise array.

    for i in range (0,len(Ti)): # Loop through all temperatures.
        temp = Ti[i] + 273.15   # Change from degrees Celsius to Kelvin
        # IAPWS package has function returning saturation pressure for given temperature in Kelvin
        Psat = (iapws97._PSat_T(temp))
        # Append pressure to array, pressure values returned in MPa
        Parray = np.append(Parray,Psat)

    # Calculate the pressure drop needed from original pressures in distribution, for boiling to occur.
    PdropNeeded = []    # Initialise array.

    for j in range (0,len(Pi)): # Loop through all pressures.
        #Pdrop needed is the difference between pressure in distribution at saturation pressure at the same temperature. 
        Pdn = Pi[j]*1.e-6 - Parray[j]
        # Append pressure to array, pressure values returned in MPa
        PdropNeeded = np.append(PdropNeeded,Pdn)

    # Find where boiling occured throughout distribution due to volcano collapse/pressure drop. 
    boilArray = []  # Initialise array.

    for k in range (0, len(PdropNeeded)):   # Loop through all pressures.
    # Boil array shows 1 if water boiled at this position and 0 if no boiling occured.
        # Water boils if Pressure drop from collapse is greater than pressure drop needed.
        if (PdropArray[k]>PdropNeeded[k]):
            boilArray = np.append(boilArray,1)
        else:
            boilArray = np.append(boilArray,0)

    # Get gold concentrations throughout distribution.
    goldSolgm3 = [] # Initialise array.

    for m in range (0,len(Ti)): # Loop through all temperatures.
        # Include if else statement for plot of contours only in boiling region.
        #if (boilArray[m] == 1):
        goldAmount = averagePorosity * exp_Fit(Ti[m], *popt)
        goldSolgm3 = np.append(goldSolgm3,goldAmount)
        #else:
        #    goldSolgm3 = np.append(goldSolgm3,0)

    DepthBoiled = []    # Initialise array.
    for i in range (0, len(rangex)):    # Loop through horizontal distances.  
        DepthBoiled = np.append(DepthBoiled,0)  # Initially append 0 for starting depth.
        for j in range(0,len(yi),61):   # Loop through depths in z direction for same horizontal distance.
            # zpoint is the current depth.
            zpoint = zi[j]
            
            if (boilArray[j+i]==1):
                # If boiling occured at this point then we find volume boiled.
                DepthBoiled = np.append(DepthBoiled,zpoint)
                # Find temperature at this point so we can get the gold concentration here.
                depthT = findTemp(rangex[i],zpoint)
                goldConc = exp_Fit(depthT,*popt)
                # Change in height is calculated as depth boiled array will always just have a start and end value.
                delZ = DepthBoiled[1]-DepthBoiled[0]
                if (i == 0):    # For first column, there is no previous rangex value so volume calculated separately.
                    volume = math.pi * (-1) *delZ * ((rangex[i])**2)    # Radial formula for cylinder.
                    minVolume = volume * minPorosity   # in m^3
                    maxVolume = volume * maxPorosity   # in m^3
                    averageVolume = volume * averagePorosity   # in m^3     
                else:   # For all other columns besides first.
                    volume = math.pi * (-1) *delZ * ((rangex[i])**2-(rangex[i-1])**2)   # Radial formula with difference in x direction also.
                    minVolume = volume * minPorosity   # in m^3
                    maxVolume = volume * maxPorosity   # in m^3
                    averageVolume = volume * averagePorosity   # in m^3 
                # Add gold amounts to totals, as volume * concentration at temperature of bottom point.
                minGold += minVolume * goldConc   # in grams
                maxGold += maxVolume * goldConc   # in grams
                averageGold += averageVolume * goldConc   # in grams

                DepthBoiled = []    # Reset depth array
                DepthBoiled = np.append(DepthBoiled,zpoint) # Append end depth as starting depth for next calculation.

            else:   # For when boiling does not occur
                DepthBoiled = []    # Reset depth array
                if (zpoint != -3800):
                    # If we have not reached the bottom of the distribution then append end depth as starting depth for next calculation.
                    DepthBoiled = np.append(DepthBoiled,zpoint) 

    # If desired, array of gold results.
    goldDeposit = [minGold, averageGold, maxGold]

    # Below is code for plotting distributions.
    '''
    #PLOTS OF PRESSURES AND DISTRIBUTIONS
    # Set up figure axes.
    ax = plt.axes([0.1,0.1,0.8,0.8])
    
    # all arrays need to be reshaped to plot correctly alongside each other.
    yi2 = yi.reshape(shp)
    # Filled contour plot of saturation pressure throughout distribution.
    #plt.contourf(yi.reshape(shp),zi.reshape(shp),Parray.reshape(shp),cmap="jet")

    # Filled contour plot of original pressure throughout distribution.
    #plt.contourf(yi.reshape(shp),zi.reshape(shp),Pi.reshape(shp),cmap="jet")

    # Filled contour plot of pressure drop required throughout distribution.
    #plt.contourf(yi.reshape(shp),zi.reshape(shp),PdropNeeded.reshape(shp),cmap="jet")

    # Contour line of boiling array showing where boiling occured throughout distribution.
    CM1 = plt.contour(yi.reshape(shp),zi.reshape(shp),boilArray.reshape(shp),colors='white',linewidths=3, linestyles='solid', levels = [0,1])   # levels = [0,0.5,1]
    
    # Filled contour plot of pressure drop from volcano collapse throughout distribution.
    #plt.contourf(yi.reshape(shp),zi.reshape(shp),PdropArray.reshape(shp),cmap="jet",levels = np.linspace(0,int(hmax/100),int(hmax/100)+1))
    
    # Contour line of gold concentrations at different temperatures throughout distribution.
    CM2 = plt.contour(yi.reshape(shp),zi.reshape(shp),goldSolgm3.reshape(shp),cmap="gray", levels = [1,2,3,4,5,6,7,8,9], hatches=['-', '/', '\\', '//'])    # Hatches if contourf
    # Label contours, manual = true if you want to manually add labels in positions of your choice
    plt.clabel(CM2, CM2.levels, inline=True, manual=True, fmt = '%.0fg/m\u00b3', fontsize=9)

    # Filled contour plot of temperatures throughout distribution.
    plt.contourf(yi.reshape(shp),zi.reshape(shp),Ti.reshape(shp),cmap="jet")
    
    # Plotting the volcano lines above distribution.
    #x = np.linspace(0,xmax,1000)
    x = np.linspace(0,xmax+50,1000)
    x2 = np.linspace(0,3*xmax/5+50,1000)
    x3 = np.linspace(3*xmax/5+50,xmax+50)
    y3 = -25 + 0*x3
    y = (-(hmax/xmax)*x+hmax)
    y2 = (-(hmax/xmax)*x2+3*hmax/5)
    plt.plot(x,y,'-.k',label='Volcano')
    plt.plot(x2,y2,':k',label='Collapsed Volcano') 
    plt.plot(x3,y3,':k')
    plt.legend()

    # Fill above distribution so that volcano lines look nicer.  
    plt.fill_between([0,3000],[0,0],[hmax+100,hmax+100], color=[0.9,0.9,0.9])
    plt.fill_between(x2,0,y2, color=[0.6,0.6,0.6])

    x1,x2,y1,y2 = plt.axis()    # Find current axis limits if needed.
    plt.axis((0,2000,-2500,hmax+100))   # Change to prefereed axis limits.

    # Show gold amount on plot
    string = "Avg Gold = {0:.2f}t".format(averageGold/1000000)
    ax.text(.70, .95, string, position = (0.58,0.81), fontstyle='italic', transform=ax.transAxes, size=11, color = 'k')
    
    # Show colour bar with label, depending on which filled contours are being plotted.
    cbar = plt.colorbar()
    #cbar.set_label('Temperature (\xb0C)')
    #cbar.set_label('Pressure Drop (MPa)')
    
    # Labels of axis and title.
    #plt.title('Pressure Drop due to Volcano Collapse')
    #plt.ylabel('Z (m)')
    #plt.xlabel('Y (m)')
    
    # Too change background colour of plots from white include.
    plt.gcf().patch.set_facecolor([255/255.,153/255.,0])
    plt.savefig('flack.png', facecolor = plt.gcf().get_facecolor(), transparent = True)
    '''

    return averageGold

# Uncomment print to run function.
print(get_Gold(1400, 700))

# Code below plots graphs of size of volcano vs amount of gold produced 
'''
# Initialise arrays
pmsL = []
pmsH = []
pms = []
size = np.linspace(200,3000,100)
avgH = 700
for sz in size: # Loop through all sizes
    pmL = get_Gold(sz,avgH)/1000000 # For change in length effects
    pmsL.append(pmL)
    pmH = get_Gold(avgH,sz)/1000000 # For change in height effects
    pmsH.append(pmH)
    pm = get_Gold(sz,sz)/1000000    # For changes in both height and length simultaneously
    pms.append(pm)     

plt.plot(size,pms,'m')
#plt.plot(size,pmsH,'r',label='Height')
#plt.title('Effect of Dimension of Volcano on\nAmount of Gold Produced',size=16)
#plt.xlabel('Dimension (m)',size=16)
#plt.ylabel('Gold (Tonnes)',size=16)
plt.legend()
#plt.show()
plt.gcf().patch.set_facecolor([255/255.,153/255.,0])
plt.savefig('sizeVSgold2w.png')#, facecolor = plt.gcf().get_facecolor(), transparent = True) 
'''    

