from numpy import *
import numpy as np
import matplotlib.pyplot as plt


def makeAry(size, waveWidth, r1, typ):
    xval = np.linspace(-waveWidth, waveWidth, size)
    yval = np.linspace(-waveWidth, waveWidth, size)
    X, Y = meshgrid(xval, yval)
    Z = arange(len(X)*len(Y)).reshape(len(X), len(Y))
    for t in range(len(X)):
        for r in range(len(Y)):
            if typ == 0:
                if (abs(X[t][r])>waveWidth*1 or abs(Y[t][r])>waveWidth*1) or (abs(X[t][r])<r1 and abs(Y[t][r])<r1):
                    Z[t][r] = 0
                else:
                    Z[t][r] = 1
            if typ == 1:
                if sqrt(X[t][r]**2+Y[t][r]**2)>waveWidth*1 or sqrt(X[t][r]**2+Y[t][r]**2)<r1*1.5:
                    Z[t][r] = 0
                else:
                    Z[t][r] = 1
    Z = np.append(Z, waveWidth)
    return Z


def Condenser(ary, fc, wavLen):
    k = 2*pi/wavLen
    waveWidth = ary[ary.size-1]
    ary = np.delete(ary, ary.size-1)
    ary = np.reshape(ary, (sqrt(ary.size), sqrt(ary.size)))
    xval = np.linspace(-waveWidth, waveWidth, sqrt(ary.size))
    yval = np.linspace(-waveWidth, waveWidth, sqrt(ary.size))
    X, Y = meshgrid(xval, yval)
    ary2 = exp(1j*k*(X**2+Y**2)/(2*fc))*ary
    ary2 = np.append(ary2, waveWidth)
    return ary2


def Propagate(inputwav, z, wavlen, areaV, dx):
    waveWidth = inputwav[inputwav.size-1]
    inputwav = np.delete(inputwav, inputwav.size-1)
    inputwav = np.reshape(inputwav, (sqrt(inputwav.size), sqrt(inputwav.size)))
    k = 2*pi/wavlen
    def Int(x, y, psi, xp, yp):
        r = sqrt(z**2+(x-xp)**2+(y-yp)**2)
        return (exp(-1j*k*r)/(r**2))*psi
        #*exp(1j*k*(x**2+y**2)/(2*fc))
    def makeInput(waval):
        xval = np.linspace(-waveWidth, waveWidth, sqrt(waval.size))
        yval = np.linspace(-waveWidth, waveWidth, sqrt(waval.size))
        X, Y = meshgrid(xval, yval)
        return X, Y
    def PsiPrimeSqr(aryx, aryy, aryP, xp, yp):
        psipr = Int(aryx, aryy, aryP, xp, yp)
        aryy = []
        aryx = []
        PsiTot = sum(psipr)
        return abs((PsiTot*z)/(1j))
    aryLx, aryLy = makeInput(inputwav)
    xp = np.linspace(-waveWidth*areaV, waveWidth*areaV, dx)
    yp = np.linspace(-waveWidth*areaV, waveWidth*areaV, dx)
    Xp, Yp = meshgrid(xp, yp)
    zs = np.array([PsiPrimeSqr(aryLx, aryLy, inputwav, xp, yp) for xp, yp in zip(np.ravel(Xp), np.ravel(Yp))])
    maxZ = np.amax(zs)
    zs = zs/maxZ
    Zp = zs.reshape(Xp.shape)
    return Xp, Yp, Zp


fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
waveWidth = .000195
r1 = .000115
waveLen = 1.5e-10
ary = makeAry(400, waveWidth, r1, 1)
'''with open('ring200_00009_0D2.txt') as f:
    for line in f:
        ary = np.array([map(float, line.split())])
    f.close()'''
aryLensed = Condenser(ary, 5, waveLen)
X, Y, Z = Propagate(aryLensed, 5, waveLen, .5, 50)
ax.plot_surface(X, Y, Z, cmap=plt.cm.jet, rstride=1, cstride=1, linewidth=0)
plt.show()
