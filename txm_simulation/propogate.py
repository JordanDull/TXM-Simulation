from numpy import *
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D


def makeAry(size, waveWidth, r1):
    """
    Makes an array that represents the plane wave. Simulates the aperature. 

    Parameters
    ----------
    size : int
        Number of subdivisions of array

    waveWidth : float
        Width of plane wave and outer radius of ring.

    r1 : float
        Inner radius of ring

    Returns
    -------
    X : float array
        Array of values on the x axis. Together with y values, they correspond to a value in Z
    Y : float array
        Array of values on the y axis. Together with x values, they correspond to a value in Z
    Z : int array
        Array of amplitudes. 1 means light is present. 0 means no light.
    """
    xval = np.linspace(-waveWidth, waveWidth, size)
    yval = np.linspace(-waveWidth, waveWidth, size)
    X, Y = meshgrid(xval, yval)
    Rsq = X**2 + Y**2
    Z = (Rsq <= (waveWidth)**2) * (Rsq >= (r1)**2)
    return X, Y, Z


def Condenser(X, Y, Z, fc, wavLen):
    """
    Simulates light passing through an infinitesimally thin lens with a focal length fc. 
    
    Parameters 
    ----------
    X : float array
        Array of values on the x axis. Together with y values, they correspond to a value in Z
        
    Y : float array
        Array of values on the y axis. Together with x values, they correspond to a value in Z
        
    Z : complex array
        Array of phases of light. 
        
    fc : float
        Focal length of lens
    
    waveLen : float
        Wavelength of light
    
    Returns
    -------
     : complex array
    Array of complex numbers simualting the phases of light at certain positions. 
    """
    '''rand=arange(X.size,dtype=complex)
    for t in range(rand.size):
        rand[t]=.0001*1j*random.random()
    rand=np.reshape(rand,(sqrt(rand.size),sqrt(rand.size)))'''
    k = 2*pi/wavLen
    return exp(1j*k*(X**2+Y**2)/(2*fc))*Z  #-rand*


def Propagate(X,Y,Z, dist, wavlen, areaV, dx):
    """
    Simulates the propagation of light given a certain distance and wavelength.
    
    Parameters
    ----------
    X : float array
        Array of values on the x axis. Together with y values, they correspond to a value in Z
        
    Y : float array
        Array of values on the y axis. Together with x values, they correspond to a value in Z
        
    Z : complex array
        Array of phases of light. 
        
    dist : float
        Distance for light to travel
        
    waveLen : float
        Wavelength of light
        
    areaV : float
        Multiplier of original width of input array. Used to form width of 
        new array. Let it equal 1 for constant width of viewing area.
    
    dx : int
        Number of subdivisions in of new array.
        
    Returns
    -------
    Xp : float array
        Array of values on the x axis. Together with y values, they correspond to a value in Zp
        
    Yp : float array
        Array of values on the y axis. Together with x values, they correspond to a value in Zp
        
    Zp : complex array
        Array of phases of light. 
    """

    def PsiPrimeSqr(aryx, aryy, aryP, xp, yp):
        """
        Calculates the amplitude of the light at a point on the new plane
        as a result of all of the values from the input array.
        
        Parameters
        ----------
        aryx : float array
            Array of values on the input x axis. Together with y values, they correspond to a value in Zp
        
        aryy : float array
            Array of values on the input y axis. Together with x values, they correspond to a value in Zp
        
        aryP : complex array
            Array of phases of light.
        
        xp : float
            x position of point being calculated for new array. 
            
        yp : float
            y position of point being calucalted for new array. 
            
        Returns
        -------
         : float
            Value of amplitude at point indicated by xp,yp. 
        """
        r = sqrt(dist**2+(aryx-xp)**2+(aryy-yp)**2)
        psipr = (exp(-1j*k*r)/(r**2))*aryP
        aryy = []
        aryx = []
        PsiTot = sum(psipr)
        return abs((PsiTot*dist)/(1j))
        
    waveWidth=abs(X[0][0])
    k = 2*pi/wavlen    
    xp = np.linspace(-waveWidth*areaV, waveWidth*areaV, dx)
    yp = np.linspace(-waveWidth*areaV, waveWidth*areaV, dx)
    Xp, Yp = meshgrid(xp, yp)
    zs = np.array([PsiPrimeSqr(X, Y, Z, xp_, yp_)
                   for xp_, yp_ in zip(np.ravel(Xp), np.ravel(Yp))])
    maxZ = np.amax(zs)
    zs = zs/maxZ
    Zp = zs.reshape(Xp.shape)
    return Xp,Yp,Zp


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
waveWidth = 195e-6
r1 =115e-6
waveLen = 1.5e-10
X,Y,Z = makeAry(100, waveWidth, r1)
'''with open('ring200_00009_0D2.txt') as f:
    for line in f:
        ary = np.array([map(float, line.split())])
    f.close()'''
ZLensed = Condenser(X, Y, Z, .05, waveLen)
Xp, Yp, Zp = Propagate(X, Y, ZLensed, .05, waveLen, .0005, 60)
ax.plot_surface(Xp, Yp, Zp, cmap=plt.cm.jet, rstride=1, cstride=1, linewidth=0)
plt.show()
