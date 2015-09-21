__author__ = 'juncheng'

import commons
import deflection
import numpy as np
import math

imgRes = 0.05
srcRes = 0.15
srcDimX=50
srcDimY=50

def testWholeGridMapping(const):
    imgXm, imgYm = commons.getGrid(-50, 50, 2, -50,  50,2, 0.05)
    deflXm = []
    deflYm = []

    for i in range(srcDimX):
        deflx, defly, _, _ = deflection.getDeflection(imgXm[i], imgYm[i],'SIE',const)
        deflXm.append(deflx)
        deflYm.append(defly)
    new_deflXm = np.asarray(deflXm)
    new_deflYm = np.asarray(deflYm)
    srcXm = imgXm - new_deflXm
    srcYm = imgYm - new_deflYm

    commons.plotGrid(srcXm, srcYm)


def test_generateEllGaussianImage(fitsName):
    n = 100
    data = np.zeros((n,n))
    s0 = 0
    s1 = 1000
    phi = np.pi/4    # 45 degree
    sigma1 = 4.95
    sigma2 = 5.00
    A = (np.cos(phi)/sigma1)**2  + (np.sin(phi)/sigma2)**2
    B = (np.sin(phi)/sigma1)**2  + (np.cos(phi)/sigma2)**2
    C = 2*np.sin(phi)*np.cos(phi)*(1/(sigma1**2)-1/(sigma2**2))
    print C

    x0 = n/2.0
    y0 = n/2.0
    for i in range(n):
        for j in range(n):
            data[i][j] = s0+s1*np.exp(-0.5*(A*(i-x0)**2 + B*(j-y0)**2 +C*(i-x0)*(j-y0)))

    commons.writeFitsImage(data, fitsName )

    #results = astro.fit_gauss_elliptical([0,0], data)
#    print results

def main():
    #test_generateEllGaussianImage("ellGaussianTest_0.01.fits")
    print commons.getTriWeight([0,0], [0, 3], [4, 0], [2,1])
    return 0

    #testWholeGridMapping(const)

if __name__ =="__main__":
    main()