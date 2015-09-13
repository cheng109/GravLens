__author__ = 'juncheng'

import commons
import deflection
import numpy as np

imgRes = 0.05
srcRes = 0.15
srcDimX=50
srcDimY=50

def testWholeGridMapping():
    imgXm, imgYm = commons.getGrid(-50, 50, 2, -50,  50,2, 0.05)
    deflXm = []
    deflYm = []

    for i in range(srcDimX):
        deflx, defly, _, _ = deflection.getDeflection(imgXm[i], imgYm[i],'SPEMD')
        deflXm.append(deflx)
        deflYm.append(defly)
    new_deflXm = np.asarray(deflXm)
    new_deflYm = np.asarray(deflYm)
    srcXm = imgXm - new_deflXm
    srcYm = imgYm - new_deflYm

    commons.plotGrid(srcXm, srcYm)


def main():
    testWholeGridMapping()

if __name__ =="__main__":
    main()