from numpy import *
import scipy
import numpy as np
import matplotlib.pyplot as plt 

def makeAry(size,waveWidth=125e-6):   
    """
    Makes an array that represents the plane wave. Simulates the aperature. 

    Parameters
    ----------
    size : int
        Number of subdivisions of array

    waveWidth : float
        Width of plane wave and outer radius of ring. Defualt is 125e-6 meters

    Returns
    -------
    X : float array
        Array of values on the x axis. Each value corresponds to a value in Z
    Z : int array
        Array of amplitudes. 1 means light is present. 0 means no light.
    """   
    X=np.linspace(-waveWidth,waveWidth,size)
    Z=(abs(X)<=waveWidth)
    return X,Z

def PropagatePoint(X,Z,dist,waveLen,areaV,dx):
    """
    Propagates a single point of light given a distance and wavelength.
    
    Parameters
    ----------
    X : float
        Location of point to be propagated
        
    Z : int 
        Amplitude of single point of light 
        
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
        Array of values on the x axis. Each corresponds to a value in Zp
                
    Zp : complex array
        Array of phases of light.
    """
    k=2*pi/waveLen
    Xp=np.linspace(-waveWidth*areaV,waveWidth*areaV,dx) 
    r=sqrt(dist**2+(X-Xp)**2)
    Zp=(exp(-1j*k*r)/r)
    return Xp,Zp
    
def Psi(X,Xp,Z,k,dist):
    """
    Calculates the amplitude of the light at a point on the new plane
    as a result of all of the values from the input array.
    
    Parameters
    ----------
    X : float array
        Array of values on the input x axis. Each value corresponds to a value in Zp
    
    Xp : float
        x position of point being calculated for new array. 
        
    Returns
    -------
     : complex
        Value of probability amplitude at point indicated by xp through summation. 
    """
    r=sqrt(dist**2+(X-Xp)**2)
    PsiP=(exp(-1j*k*r)/r)*Z
    return sum(PsiP)
    
def Propagate(X,Z,dist,waveLen,areaV,dx,shift):
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
    
    shift : float
        Shifts the detection area by distance indicated. 
        
    Returns
    -------
    Xp : float array
        Array of values on the x axis. Each corresponds to a value in Zp
                
    Zp : complex array
        Array of phases of light. 
    """
    k=2*pi/waveLen
    Xp=np.linspace(-waveWidth*areaV+shift,waveWidth*areaV+shift,dx)
    Zp=np.array([Psi(X,xp_,Z,k,dist) for xp_ in Xp])
    return Xp,Zp
    
def Condenser(X,Z,f,wavLen,apr,shift):
    """
    Simulates light passing through an infinitesimally thin lens with a focal length f. 
    
    Parameters 
    ----------
    X : float array
        Array of values on the x axis. Each value corresponds to a value in Z
        
    Z : complex array
        Array of phases of light. 
        
    f : float
        Focal length of lens
    
    waveLen : float
        Wavelength of light
    
    apr : int
        Indicates which aperature to use. 1 for condenser or 0 for objective
        
    shift : float
        Shifts the objective lens by distance indicated
    
    Returns
    -------
     : complex array
        Array of complex numbers simualting the phases of light at certain positions. 
    """ 
    k=2*pi/wavLen
    if apr==1:
        Z=(Z*(X<=waveWidth)*(X>=r1))+(1*Z*(-X<=waveWidth)*(-X>=r1)) #Z*(abs(X)<=waveWidth)*(abs(X)>=r1)
        return exp(1j*k*(X**2)/(2*f))*Z
    if apr==0:
        Z=Z*(X<=50e-6+shift)*(X>=-50e-6+shift)   #Z*(X<=50e-6)*(X>=-50e-6)
        return exp(1j*k*((X-shift)**2)/(2*f))*Z          #  !!!!!!!!!!!!!!!!!!!!!!!!!

def Sample(X,Z,index,waveLen):
    """
    Simulates light passing through a sample of a given index of refraction, 1 micron thick. 
    
    Parameters 
    ----------
    X : float array
        Array of values on the x axis. Each value corresponds to a value in Z
        
    Z : complex array
        Array of phases of light. 
        
    index : complex
        Index of refraction
    
    waveLen : float
        Wavelength of light
    
    Returns
    -------
    ary2 : complex array
        Array of complex numbers simualting the phases of light at certain positions. 
    """
    k=2*pi/waveLen 
    sample=(abs(X)<=5e-6)*1e-6
    ary2=Z*exp(-1j*k*index*(sample)) 
    return ary2

#Constants    
waveWidth=195e-6
r1=115e-6 
waveLen=1.5e-10
dist=1.9
fc=.1
fo=.0194
n=1-(4.77e-5)-1j*(4.95e-6)
aprCon=1
aprObs=0

#Summation Approach
X,Z=makeAry(700) #700
size=2500 #500
ans=arange(X.size*size).reshape(X.size,size)
ansSum=[0]*size
for i in range(X.size):
    Xp,Zp=PropagatePoint(X[i],Z[i],dist,waveLen,1,3000)
    ZL=Condenser(Xp,Zp,fc,waveLen,aprCon,0)
    Xpp,Zpp=Propagate(Xp,ZL,.105556,waveLen,.1,2500,0)
    ZSam=Sample(Xpp,Zpp, n, waveLen)
    Xppp,Zppp=Propagate(Xpp,ZSam,.01959,waveLen,.5,size,0) #.01959   #2500
    #ZpppL=Condenser(Xppp,Zppp/100,fo,waveLen,aprObs,0)
    #X4p,Z4p=Propagate(Xppp,ZpppL,2,waveLen,6,size,0)
    ans[i]=abs(Zppp)
    ansSum=ansSum+(ans[i]/1000)**2
    
'''#Single Point Propagation
X,Z=-.00009,1  #.00009
Xp,Zp=PropagatePoint(X,Z,dist,waveLen,1,3000)
ZL=Condenser(Xp,Zp,fc,waveLen,aprCon)
Xpp,Zpp=Propagate(Xp,ZL,.105556,waveLen,.1,2500)
ZSam=Sample(Xpp,Zpp, n, waveLen)
Xppp,Zppp=Propagate(Xpp,ZSam,.01959,waveLen,.5,2500)
ZpppL=Condenser(Xppp,Zppp,fo,waveLen,aprObs)
X4p,Z4p=Propagate(Xppp,ZpppL,2,waveLen,6,1000)'''

a=open('sum2_obs_sam.txt','w')
for t in range(ansSum.size):
    a.write(str(ansSum[t])+' ')
a.close()

fig=plt.figure()
#plt.plot(X4p,abs(Z4p))
plt.plot(Xppp,ansSum)#,X4p,ans.T)
plt.ylabel('Light Intensity')
plt.xlabel('Position (m)')

plt.show()
