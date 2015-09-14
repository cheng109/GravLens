__author__ = 'juncheng'

import commons
import deflection
import numpy as np

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


def main():
    dirName= "test_images/"
    imgData = commons.readFitsImage(dirName + "model_img.fits")
    srcData = commons.readFitsImage(dirName + "model_src.fits")
    varData = commons.readFitsImage(dirName+ "var.fits.gz")
    psfData = commons.readFitsImage(dirName+ "psf.fits.gz")
    const = commons.Constants(srcData.shape,imgData.shape,0.05, 0.05)

    #testWholeGridMapping(const)

if __name__ =="__main__":
    main()