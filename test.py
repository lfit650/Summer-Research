# Function to find temperature at a point y,z underground in a
# Geothermal System.
#import import_Temperature
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

# read temperature and pressure data 
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

# dy = zi[0]
# for i in range (1,len(zi)):
#     if (zi[i]<dy):
#         dy = zi[i]
# print(dy)
#print(zi)
#print(yi)
##print(zi)
#print(Ti)
#print(Pi)
#for i in range(61):
#    print(zi[i])
#print(len(Ti))
#print(len(zi))

'''
# TEMPERATURE DISTRIBUTION

yi2 = yi.reshape(shp)
plt.contourf(yi.reshape(shp),zi.reshape(shp),Ti.reshape(shp),cmap="jet")
plt.title('Temperature Distribution in Geothermal System')
plt.xlabel('Y (m)')
plt.ylabel('Z (m)')
x1,x2,y1,y2 = plt.axis()
plt.axis((0,2000,y1,y2))
# sortedT = np.sort(Ti)
# CS = plt.imshow(sortedT, origin='lower', cmap=cm.jet, extent = [yi.min(),yi.max(),zi.min(),zi.max()], aspect='auto')


cbar = plt.colorbar()
cbar.set_label('Temperature (\xb0C)')
plt.savefig('FinalTempDist2000.png')
'''


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

def pressureDrop(y,xmax,hmax):
	# hmax = 				# max height of volcano
    #xmax = 0.5*hmax				# max radius, distance from centre to edge of volcano. 
    delh = -(hmax/xmax) * y + hmax
    Pdrop = Pressure(delh)
    return Pdrop  

def Phase(P,T):
	#Tc = 647 # Kelvin
	#Pc = 22.064 # MPa
	P = P/1000000 # Change to MPa
	#return _utils.getphase(Tc = Tc, Pc = Pc, T = T, P = P)
	ph = IAPWS97(T=T, P=P)
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


# Gives solubility data due to change in temperature only at different pH levels. Pressure not measured/constant
# For graph c pH at 20 degrees C of 7.45.
xTemp = np.array([165, 200, 213.333333, 230, 242, 252, 268.5, 279.8, 290, 296])
yConc = np.array([5, 10, 15, 20, 25, 30, 40, 50, 60, 67])
# Check because gives negative values at low temperatures.
    
# For graph a, pH at 20 degrees C of 9.55.
# xTemp = np.array([164.9, 190, 240, 260, 268, 280, 284.5, 290, 298.3, 300])
# yConc = np.array([2.5, 3, 5, 8.5, 10, 13.5, 15, 17, 20, 21.25])

def exp_Fit(x, a, b, c):
    return a*np.exp(b*x)+c

popt, pcov = curve_fit(exp_Fit, xTemp, yConc, p0=(1, 1e-6, 1))
'''
xxTemp = np.linspace(0, 320, 2501)
yyConc = exp_Fit(xxTemp, *popt)

plt.plot(xTemp,yConc,'o', xxTemp, yyConc)
plt.title('Gold Solubility Changes due to Temperature')
plt.xlabel('Temp. (\xb0C)')
plt.ylabel('Au Conc. (ppm)')
plt.savefig('expFit for c.png')
'''
# Call exp_Fit(temperatureValue) ?? exp_Fit(temperatureValue, *popt) to get concentration.


def findPressure(y,z,xmax,hmax):
    # # hmax = max height of volcano

    # # Change xmax for different type of volcanoes (angles).
    # xmax = hmax # max radius, distance from centre to edge of volcano. 45 degrees if = hmax
    # delh = -(hmax/xmax) * y + hmax
    # Pi = Pressure(delh-z)
    # return Pi
    #xmax = 0.5*hmax # max radius, distance from centre to edge of volcano. 45 degrees if = hmax
    delh = -(hmax/xmax) * y + hmax
    i = 0
    while (zi[i] != z):
        i+=1
    j = i
    while ((yi[j]) != y and j<i+61):
        j+=1
    P = Pi[j]
    P = P + Pressure(delh)
    return P

def testPhase(Pi, Ti):
    if (Phase(Pi,Ti) != "Vapour"):
        return False
    return True
#print(Phase(356744.375, 402.039))
def get_Gold(xmax,hmax):

    minPorosity = 0.045
    maxPorosity = 0.0966666667
    averagePorosity = 0.07083333333333333    
    minGold = 0
    maxGold = 0
    averageGold = 0
    # hmax = 
    #x1 = hmax
    # n = 1000
    # rangex = np.linspace(0,x1,num = n, endpoint = False)
    # Pdrop = [0]*n
    # for i in range(len(rangex)):
    #     np.append(Pdrop,pressureDrop(rangex[i],hmax)
    #     #print(Pdrop[i])

    j = 0 
    while (xmax >= yi[j]):   # Find values we can use from temperature distribution
        j+=1
    # Initialise arrays    
    rangex = yi[0:j]   
    #print(rangex) 
    
    # CHANGE THIS FOR FLANK MODEL FLAT TO HALFWAY OR JUST BEFORE --> PARRALEL TO CENTRE --> VERTICALAT CENTRE

    Pdrop = [0]*2501
    for i in range(0, len(rangex)): # Loop through sections in x direction
            j=i
            Pdrop[j] = (1.e-6)*pressureDrop(rangex[i],xmax,hmax)
            while (j+61<len(yi)):
                j+=61
                Pdrop[j] = (1.e-6)*pressureDrop(rangex[i],xmax,hmax)
    PdropArray = np.asarray(Pdrop)
    #print(Pdrop)
    # Pressure required
    Parray = []
    for i in range (0,len(Ti)):
        temp = Ti[i] + 273.15
        Psat = (iapws97._PSat_T(temp))
        Parray = np.append(Parray,Psat)

    #print(Parray)
    PdropNeeded = []
    for j in range (0,len(Pi)):
        Pdn = Pi[j]*1.e-6 - Parray[j]
        PdropNeeded = np.append(PdropNeeded,Pdn)

    boilArray = []
    for k in range (0, len(PdropNeeded)):
        if (PdropArray[k]>PdropNeeded[k]):
            boilArray = np.append(boilArray,1)
        else:
            boilArray = np.append(boilArray,0)

    # for i in range (0,len(boilArray)):
    #     print(boilArray[i])
    goldSolgm3 = []
    for m in range (0,len(Ti)):
        # Include if else statement for contours only in boiling region.
        #if (boilArray[m] == 1):
        goldAmount = averagePorosity * exp_Fit(Ti[m], *popt)
        goldSolgm3 = np.append(goldSolgm3,goldAmount)
        #else:
        #    goldSolgm3 = np.append(goldSolgm3,0)

    #print(goldSolgm3)
    amountGold = [0]*2501
    region = []
    for n in range (0,len(Ti)):
        region = np.append(region,iapws97._Bound_TP(Ti[n]+273.15,(1.e-6)*Pi[n]))
    
    # To use gold solubility changes need to calculate all cylinders at every depth and find concentration for bottom of the cylinder temperature.
    '''
    # NO CHANGE IN SOLUBILITY

    DepthBoiled = []
    totalVolume = 0
    for i in range (0, len(rangex)):
        DepthBoiled = np.append(DepthBoiled,0)
        x = True
        for j in range(0,len(yi),61):   #Loop through depths in z direction
            zpoint = zi[j]
            if ((boilArray[j+i]==1) and (x == False)):
                x=True
                DepthBoiled = np.append(DepthBoiled,zpoint)
            if (x==True):
                if (boilArray[j+i]==0):
                    zpoint = zi[j-61]
                    DepthBoiled = np.append(DepthBoiled,zpoint)
                    delZ = DepthBoiled[1]-DepthBoiled[0]
                    if (i == 0):
                        volume = math.pi * (-1) *delZ * ((rangex[i])**2)
                    else:    
                        volume = math.pi * (-1) *delZ * ((rangex[i])**2-(rangex[i-1])**2)
                    totalVolume += volume
                    volume = 0
                    DepthBoiled = []
                    x=False
    
    minVolume = totalVolume * minPorosity   # in m^3
    maxVolume = totalVolume * maxPorosity   # in m^3
    averageVolume = totalVolume * averagePorosity   # in m^3 
    
    minGoldConc = 4.78358   # in ppm
    maxGoldConc = 16.33333  # in ppm
    averageGoldConc = 10.55845833333333 # in ppm        # Subject to change at different depths, graphs in parameter research.

    minGold = minVolume * minGoldConc   # in grams
    maxGold = maxVolume * maxGoldConc   # in grams
    averageGold = averageVolume * averageGoldConc   # in grams   
    #print(DepthBoiled)
    #print(zi)
    '''
    # CHANGE IN SOLUBILITY INCLUDED
    #zPointArray = []
    DepthBoiled = []
    for i in range (0, len(rangex)):
        DepthBoiled = np.append(DepthBoiled,0)
        for j in range(0,len(yi),61):   #Loop through depths in z direction
            zpoint = zi[j]
            #print(j+i)
            #print(boilArray[j+i])
            if (boilArray[j+i]==1):
                DepthBoiled = np.append(DepthBoiled,zpoint)
                depthT = findTemp(rangex[i],zpoint)
                goldConc = exp_Fit(depthT, *popt)
                delZ = DepthBoiled[1]-DepthBoiled[0]
                if (i == 0):
                    volume = math.pi * (-1) *delZ * ((rangex[i])**2)
                    minVolume = volume * minPorosity   # in m^3
                    maxVolume = volume * maxPorosity   # in m^3
                    averageVolume = volume * averagePorosity   # in m^3     
                else:    
                    volume = math.pi * (-1) *delZ * ((rangex[i])**2-(rangex[i-1])**2)
                    minVolume = volume * minPorosity   # in m^3
                    maxVolume = volume * maxPorosity   # in m^3
                    averageVolume = volume * averagePorosity   # in m^3 

                minGold += minVolume * goldConc   # in grams
                maxGold += maxVolume * goldConc   # in grams
                averageGold += averageVolume * goldConc   # in grams
                amountGold[j+i] = averageGold

                DepthBoiled = []
                DepthBoiled = np.append(DepthBoiled,zpoint)

            else:
                DepthBoiled = []
                if (zpoint != -3800):
                    DepthBoiled = np.append(DepthBoiled,zpoint)

    goldDeposit = [minGold, averageGold, maxGold]

    amountGoldArray = np.asarray(amountGold)
    '''
    #PLOTS OF PRESSURES AND DISTRIBUTIONS
    ax = plt.axes([0.1,0.1,0.8,0.8])
    yi2 = yi.reshape(shp)
    #plt.contourf(yi.reshape(shp),zi.reshape(shp),Parray.reshape(shp),cmap="jet")
    #plt.contourf(yi.reshape(shp),zi.reshape(shp),Pi.reshape(shp),cmap="jet")
    #plt.contourf(yi.reshape(shp),zi.reshape(shp),PdropNeeded.reshape(shp),cmap="jet")
    CM1 = plt.contour(yi.reshape(shp),zi.reshape(shp),boilArray.reshape(shp),colors='white',linewidths=3, linestyles='solid', levels = [0,1])   # levels = [0,0.5,1]
    #plt.contourf(yi.reshape(shp),zi.reshape(shp),PdropArray.reshape(shp),cmap="jet",levels = np.linspace(0,int(hmax/100),int(hmax/100)+1))
    CM2 = plt.contour(yi.reshape(shp),zi.reshape(shp),goldSolgm3.reshape(shp),cmap="gray", levels = [1,2,3,4,5,6,7,8,9], hatches=['-', '/', '\\', '//'])    # Hatches if contourf
    plt.clabel(CM2, CM2.levels, inline=True, manual=True, fmt = '%.0fg/m\u00b3', fontsize=9)
    plt.contourf(yi.reshape(shp),zi.reshape(shp),Ti.reshape(shp),cmap="jet")
    #x = np.linspace(0,xmax,1000)
    x = np.linspace(0,xmax+50,1000)
    plt.plot(x,(-(hmax/xmax)*x+hmax),'-.k',label='Volcano')
    plt.legend()

    plt.fill_between([0,3000],[0,0],[hmax+100,hmax+100], color=[0.9,0.9,0.9])
    

    x1,x2,y1,y2 = plt.axis()
    plt.axis((0,2000,-2500,hmax+100))
    string = "Avg Gold = {0:.2f}t".format(averageGold/1000000)
    ax.text(.70, .95, string, position = (0.57,0.83), fontstyle='italic', transform=ax.transAxes, size=11, color = 'k')
    cbar = plt.colorbar()
    cbar.set_label('Temperature (\xb0C)')
    #cbar.set_label('Pressure Drop (MPa)')
    
    plt.title('Pressure Drop due to Volcano Collapse')
    plt.ylabel('Z (m)')
    plt.xlabel('Y (m)')
    
    #plt.savefig('Parray.png')
    #plt.savefig('Pi.png')
    #plt.savefig('PdropNeeded.png')
    #plt.savefig('boilArrayxmaX14.png')
    #plt.savefig('Pdrop7x14.png')
    plt.savefig('final.png')
    
    Tarray = [100]
    depthArray = [0]
    T0 = 373.15 #K
    p0 = 101325 #Pa
    pv = 0.8    #kg/m^3
    lw = 2256000    #J/kg
    pw = 1000   #kg/m^3
    Ts = [100]
    i = 3
    for q in range (0+i,len(zi),61):
        depthArray = np.append(depthArray,zi[q])
        ps = Pi[q]
        Ts = np.append(Ts,Ti[q])
        # http://physicstasks.eu/1702/boiling-point-of-water-at-high-pressure
        #Tw = T0*np.exp(((ps-p0)*(1/pv-1/pw))/lw)
        #Tc = Tw-273.15
        
        # http://sciencechatforum.com/viewtopic.php?f=83&t=4977&start=0
        #Temperature = -3200 / np.log(ps*0.000145038/1700) - 459.7
        #Tc = 5/9*(Temperature-32)
        try:
            Tc = iapws97._TSat_P(ps/1.e6)-273.15
        except:
            Tc = -1
        Tarray = np.append(Tarray,Tc)
    
    plt.plot(Ts,depthArray,'r:')
    #plt.plot(Ts,depthArray+700.,'r--')
    '''
    plt.plot(Ti,Parray,'b-')
    plt.title('Temperature VS Saturation Pressure',size=16)
    plt.xlabel('Temperature (\xb0C)',size=16)
    plt.ylabel('Sat. Pressure (MPa)',size=16)
    #x1,x2,y1,y2 = plt.axis()
    #plt.axis((100,x2,-2000,0))    
    plt.show()
    
    #plt.savefig('parrayVST.png')
    #return goldDeposit  #in grams
    return averageGold

    # Same as "low quality underground mine"


print(get_Gold(1400, 700))
'''   
pmsL = []
pmsH = []
pms = []
size = np.linspace(200,3000,100)
avgH = 700
for sz in size:
    pmL = get_Gold(sz,avgH)/1000000
    pmsL.append(pmL)
    pmH = get_Gold(avgH,sz)/1000000
    pmsH.append(pmH)
    pm = get_Gold(sz,sz)/1000000
    pms.append(pm)     
#for i in range (0,len(x)):
#for sz,pm in zip(size,pms):
plt.plot(size,pmsL,'b',label='Length')
plt.plot(size,pmsH,'r',label='Height')
plt.title('Effect of Dimension of Volcano on\nAmount of Gold Produced',size=16)
plt.xlabel('Dimension (m)',size=16)
plt.ylabel('Gold (Tonnes)',size=16)
plt.legend()
plt.show()
#plt.savefig('sizeVSgold.png')     
'''    