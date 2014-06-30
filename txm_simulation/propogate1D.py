from numpy import *
import scipy
import numpy as np
import matplotlib.pyplot as plt 
import random

def makeAry(size,waveWidth,r1): 
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
        Array of values on the x axis. Each value corresponds to a value in Z
    Z : int array
        Array of amplitudes. 1 means light is present. 0 means no light.
    """   
    X=np.linspace(-waveWidth,waveWidth,size)
    Z=(abs(X)<=waveWidth)*(abs(X)>=r1)
    return X,Z
    
def Condenser(X,Z,fc,wavLen):
    """
    Simulates light passing through an infinitesimally thin lens with a focal length fc. 
    
    Parameters 
    ----------
    X : float array
        Array of values on the x axis. Each value corresponds to a value in Z
        
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
    for t in range(X.size):
        rand[t]=.0001*1j*random.random()'''
    k=2*pi/wavLen
    return exp(1j*k*(X**2)/(2*fc))*Z  #, (-rand*1j*(X**2)/(2*fc))*Z multiply by rand in exp if want noise 
    
def Propagate(X,Z,dist,wavlen,areaV,dx): 
    """
    Simulates the propagation of light given a certain distance and wavelength.
    
    Parameters
    ----------
    X : float array
        Array of values on the x axis. Each value corresponds to a value in Z
        
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
    
    def PsiPrimeSqr(X,Z,xp):
        """
        Calculates the amplitude of the light at a point on the new plane
        as a result of all of the values from the input array.
        
        Parameters
        ----------
        aryx : float array
            Array of values on the input x axis. Each value corresponds to a value in Zp
        
        aryP : complex array
            Array of phases of light.
        
        xp : float
            x position of point being calculated for new array. 
            
        Returns
        -------
         : float
            Value of amplitude at point indicated by xp,yp. 
        """
        r=sqrt(dist**2+(X-xp)**2)
        psipr=(exp(-1j*k*r)/(r**2))*Z 
        PsiTot=sum(psipr)
        return abs((PsiTot*dist)/(1j))  
    k=2*pi/wavlen 
    waveWidth=abs(X[0])
    Xp=np.linspace(-waveWidth*areaV,waveWidth*areaV,dx)
    Zp=np.array([PsiPrimeSqr(X,Z,xp_) for xp_ in Xp])
    maxZ=np.amax(Zp)
    Zp=Zp/maxZ
    return Xp,Zp
def Sample(X,Z,index,waveLen):
    k=2*pi/waveLen 
    sample=(abs(X)<=.00001)*.000001
    ary2=Z*exp(1j*k*index*(sample)) 
    return abs(ary2/1j)
    
fig=plt.figure()   
'''with open('square1000_00009_0D1.txt') as f:
    for line in f:
        ary = np.array([map(float, line.split())])
    f.close()'''
waveWidth=195e-6 
r1=115e-6
waveLen=1.5e-10
n=1-4.77e-5-1j*4.95e-6
X,Z=makeAry(20000,waveWidth,r1)
ZLensed=Condenser(X,Z,.05,waveLen)
Xp,Zp=Propagate(X,ZLensed,.05,waveLen,.005,1000)
#Zpp=Sample(Xp,Zp,n,waveLen)
#y=np.append(y,abs(x[0]))
#aryLensed2=Condenser(y,2,waveLen)
#xp,yp=Propagate(aryLensed2,4,waveLen,.75,500)
#plt.subplot(211)
plt.plot(Xp,Zp)
plt.ylabel('Light Intensity')
plt.xlabel('Position')
#plt.subplot(212)
#plt.plot(X,Lens)

plt.show()
