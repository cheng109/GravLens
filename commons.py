
from astropy.io import fits
from fastell4py import fastell4py
import numpy as np

class SPEMD:
    def _init__(self):
        self.name="SPEMD"
        self.critRad=0
        self.ellipticity = 0
        self.axisRatio = 0
        self.pa = 0
        self.gamma = 0
        self.coreRad = 0
        self.centerX=0
        self.centerY=0

def readFitsImage(imageName):
    hdulist = fits.open(imageName)
    data =  hdulist[0].data
    hdulist.close()
    return data

def writeFitsImage(data, outputName):
    hdu = fits.PrimaryHDU(data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(outputName, clobber=True)

def test_fastelldefl():
    x1 = np.array([1,2,3])
    x2 = np.array([3,2,1])
    q = 0.9
    gam = 0.5
    arat = 0.9
    s = 0.01
    alpha1, alpha2 = fastell4py.fastelldefl(x1, x2, q, gam, arat, s)
    print alpha1, alpha2, 'alpha1, alpha2'
    assert alpha1[0] == 0.49802187030058409
    assert alpha2[0] == 1.604771432322617

    x1 = 1.
    x2 = 3.
    alpha1_new, alpha2_new = fastell4py.fastelldefl(x1, x2, q, gam, arat, s)
    assert alpha1_new == alpha1[0]
    assert alpha2_new == alpha2[0]


def getGrid(xStart, xEnd, xStep, yStart, yEnd, yStep, pixelSize):
    x = [pixelSize*x for x in np.arange(xStart, xEnd, xStep)]
    y = [pixelSize*y for y in np.arange(yStart,yEnd, yStep)]
    xm, ym = np.meshgrid(x, y)
    return xm, ym


def main():
    test_fastelldefl()
if __name__=='__main__':
    main()
