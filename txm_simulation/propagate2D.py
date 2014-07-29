from numpy import *
import scipy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from multiprocessing import Process, Queue

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
X,Y,Z=makeAry(20)
size=150

def do_sum(X,Y,Z):
  print len(X),len(Y[0])
  size=150
  ans=arange(X.size*size**2, dtype=complex).reshape(X.size,size**2)
  ansSum=np.array([0]*size**2).reshape(size,size)
  for i in range(len(X)):
    for j in range(len(Y[0])):
        Xp,Yp,Zp=PropagatePoint(X[i][j],Y[i][j],Z[i][j],dist,waveLen,1,1500)
        ZL=Condenser(Xp,Yp,Zp,fc,waveLen)
        Xpp,Ypp,Zpp=Propagate(Xp,Yp,ZL,.105556,waveLen,.1,size)
        ansSum=ansSum+abs(Zpp)
        Zpp=Zpp.reshape(1,size**2)
        ans[i*len(X)+j]=Zpp[0]
  return ansSum,ans

def do_prop (q,X,Y,Z):
  q.put(do_sum(X,Y,Z))
  
def main():
  q=Queue()
  p1=Process(target=do_prop, args=(q,X[:2][:20],Y[:2][:20],Z[:2][:20]))
  p2=Process(target=do_prop, args=(q,X[2:4][:20],Y[2:4][:20],Z[2:4][:20]))
  p3=Process(target=do_prop, args=(q,X[4:6][:20],Y[4:6][:20],Z[4:6][:20]))
  p4=Process(target=do_prop, args=(q,X[6:8][:20],Y[6:8][:20],Z[6:8][:20]))
  p5=Process(target=do_prop, args=(q,X[8:10][:20],Y[8:10][:20],Z[8:10][:20]))
  p6=Process(target=do_prop, args=(q,X[10:12][:20],Y[10:12][:20],Z[10:12][:20]))
  p7=Process(target=do_prop, args=(q,X[12:14][:20],Y[12:14][:20],Z[12:14][:20]))
  p8=Process(target=do_prop, args=(q,X[14:16][:20],Y[14:16][:20],Z[14:16][:20]))
  p9=Process(target=do_prop, args=(q,X[16:18][:20],Y[16:18][:20],Z[16:18][:20]))
  p10=Process(target=do_prop, args=(q,X[18:][:20],Y[18:][:20],Z[18:][:20]))
   
  p1.start()
  p2.start()
  p3.start()
  p4.start()
  p5.start()
  p6.start()
  p7.start()
  p8.start()
  p9.start()
  p10.start()
  
  r1=q.get()
  r2=q.get()
  r3=q.get()
  r4=q.get()
  r5=q.get()
  r6=q.get()
  r7=q.get()
  r8=q.get()
  r9=q.get()
  r10=q.get()
  
  ansSum=r1[0]+r2[0]+r3[0]+r4[0]+r5[0]+r6[0]+r7[0]+r8[0]+r9[0]+r10[0]
  ans=np.array([r1[1],r2[1],r3[1],r4[1],r5[1],r6[1],r7[1],r8[1],r9[1],r10[1]])
  return ansSum,ans
        
        
'''X,Y,Z=.000,.000,1
Xp,Yp,Zp=PropagatePoint(X,Y,Z,dist,waveLen,1,1500)
ZL=Condenser(Xp,Yp,Zp,fc,waveLen)
Xpp,Ypp,Zpp=Propagate(Xp,Yp,ZL,.105556,waveLen,.1,100)
#Xppp,Yppp,Zppp=Propagate(Xpp,Ypp,Zpp/100,.0253165, waveLen,1,20)'''


Xpp,Ypp,Z=makeAry(20)

if __name__=='__main__':
  Zpp,ans=main()
  
fig=plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(Xpp, Ypp, Zpp, cmap=plt.cm.jet, rstride=1, cstride=1, linewidth=0)
ax.set_xlabel('X Position (m)')
ax.set_ylabel('Y Position (m)')
ax.set_zlabel('Light Intensity')
plt.show()

a=open('sum2D_sam_nosam.txt','w')
#a.write(str(X.size)+' '+str(size**2)+' ')
for j in range(10):
  for t in range(40):
    for s in range(size**2):
      a.write(str(ans[j][t][s])+' ')
a.close()


