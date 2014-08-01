from numpy import *
import scipy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def makeAry(size,waveWidth=125e-6):   
    xval = np.linspace(-waveWidth, waveWidth, size)
    yval = np.linspace(-waveWidth, waveWidth, size)
    X, Y = meshgrid(xval, yval)
    Rsq = X**2 + Y**2
    Z = (Rsq <= (waveWidth)**2)
    return X, Y, Z

def PropagatePoint(X,Y,Z,dist,waveLen,areaV,dx):
    def Psi(X,Y,Z,Xp,Yp):
        r=sqrt(dist**2+(X-Xp)**2+(Y-Yp)**2)
        return (exp(-1j*k*r)/r)
    k=2*pi/waveLen
    xp=np.linspace(-waveWidth*areaV,waveWidth*areaV,dx)
    yp=np.linspace(-waveWidth*areaV,waveWidth*areaV,dx)
    Xp, Yp = meshgrid(xp, yp)
    zs = np.array([Psi(X, Y, Z, xp_, yp_)
                   for xp_, yp_ in zip(np.ravel(Xp), np.ravel(Yp))])
    Zp = zs.reshape(Xp.shape)
    return Xp,Yp,Zp

def Propagate(X,Y,Z,dist,waveLen,areaV,dx):
    def Psi(X,Y,Z,Xp,Yp):
        r=sqrt(dist**2+(X-Xp)**2+(Y-Yp)**2)
        PsiP=(exp(-1j*k*r)/r)*Z
        return sum(PsiP)
    k=2*pi/waveLen
    xp=np.linspace(-waveWidth*areaV,waveWidth*areaV,dx)
    yp=np.linspace(-waveWidth*areaV,waveWidth*areaV,dx)
    Xp, Yp = meshgrid(xp, yp)
    zs = np.array([Psi(X, Y, Z, xp_, yp_)
                   for xp_, yp_ in zip(np.ravel(Xp), np.ravel(Yp))])
    Zp = zs.reshape(Xp.shape)
    return Xp,Yp,Zp
    
def Condenser(X,Y,Z,fc,wavLen):
    Rsq = X**2 + Y**2
    Z = Z*(Rsq <= (waveWidth)**2) * (Rsq >= (r1)**2)
    k=2*pi/wavLen
    return exp(1j*k*(X**2+Y**2)/(2*fc))*Z

def Sample(X,Y,Z,index,waveLen):
    k=2*pi/waveLen 
    sample=(abs(X)<=5e-6)*(abs(Y)<=5e-6)*1e-6
    ary2=Z*exp(-1j*k*index*(sample)) 
    return ary2
    
waveWidth=195e-6 
r1=115e-6
waveLen=1.5e-10
dist=1.9
fc=.1
X,Y,Z=makeAry(30)
size=200
ans=arange(X.size*size**2, dtype=complex).reshape(X.size,size**2)
ansSum=np.array([0]*size**2).reshape(size,size)
for i in range(len(X)):
    for j in range(len(Y)):
        Xp,Yp,Zp=PropagatePoint(X[i][j],Y[i][j],Z[i][j],dist,waveLen,1,1500)
        ZL=Condenser(Xp,Yp,Zp,fc,waveLen)
        Xpp,Ypp,Zpp=Propagate(Xp,Yp,ZL,.105556,waveLen,.1,size)
        ansSum=ansSum+(abs(Zpp)/1000)**2
        Zpp=Zpp.reshape(1,size**2)
        ans[i*len(X)+j]=Zpp[0]
        
'''X,Y,Z=.000,.000,1
Xp,Yp,Zp=PropagatePoint(X,Y,Z,dist,waveLen,1,1500)
ZL=Condenser(Xp,Yp,Zp,fc,waveLen)
Xpp,Ypp,Zpp=Propagate(Xp,Yp,ZL,.105556,waveLen,.1,100)
#Xppp,Yppp,Zppp=Propagate(Xpp,Ypp,Zpp/100,.0253165, waveLen,1,20)'''

a=open('sum2D_sam_nosam.txt','w')
a.write(str(X.size)+' '+str(size**2)+' ')
for t in range(X.size):
    for s in range(size**2):
        a.write(str(ans[t][s])+' ')
a.close()

fig=plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(Xpp, Ypp, ansSum, cmap=plt.cm.jet, rstride=1, cstride=1, linewidth=0)
ax.set_xlabel('X Position (m)')
ax.set_ylabel('Y Position (m)')
ax.set_zlabel('Light Intensity')
plt.show()
