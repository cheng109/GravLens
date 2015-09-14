
from astropy.io import fits
from fastell4py import fastell4py
import numpy as np
import matplotlib.pyplot as plt


class Constants():
    def __init__(self, srcSize, imgSize, srcRes, imgRes):
        self.srcSize = srcSize
        self.imgSize = imgSize
        self.srcRes = srcRes
        self.imgRes = imgRes
        self.srcXCenter = srcSize[0]/2.0
        self.srcYCenter = srcSize[1]/2.0
        self.imgXCenter = imgSize[0]/2.0
        self.imgYCenter = imgSize[1]/2.0



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


def plotGrid(xm, ym):
    #xm, ym = commons.getGrid(xSize, ySize)
    plt.plot(xm, ym, '-b')
    plt.xlim([-0.5, 0.5])
    plt.plot(ym, xm, '-b')
    plt.ylim([-0.5,0.5])
    plt.show()


def applyMask(maskFileName, mappingDict):

    maskData = readFitsImage(maskFileName)
    newMappingDict = {}
    for i in range(maskData.shape[0]):
        for j in range(maskData.shape[1]):
            if maskData[i][j]>0.5:
                newMappingDict[(i,j)]=mappingDict[(i,j)]
    return newMappingDict


def plotMappingDict(mappingDict,const):
    #### mappingDict={'imageGrid': srcGrid,  'imageGrid':srcGrid, .....}
    f, (ax1, ax2) = plt.subplots(1, 2) #, sharex=True, sharey=True )

    imgPointList = mappingDict.keys()
    srcPointList = mappingDict.values()

    srcXm, srcYm = getGrid(xStart=-2, xEnd=13, xStep=1, yStart=-2, yEnd=13, yStep=1, pixelSize=1)
    ax1.plot(srcXm, srcYm, 'r-', linewidth=0.5)
    ax1.plot(srcYm, srcXm, 'r-', linewidth=0.5)


    for i in range(len(imgPointList)):
        for j in np.arange(i+1, len(imgPointList), 1):
            if (imgPointList[i][0]==imgPointList[j][0] and abs(imgPointList[i][1]-imgPointList[j][1])==1) or (imgPointList[i][1]==imgPointList[j][1] and abs(imgPointList[i][0]-imgPointList[j][0])==1):
                ax1.plot((srcPointList[i][0], srcPointList[j][0]),(srcPointList[i][1], srcPointList[j][1]) , 'b-')
                ax2.plot((imgPointList[i][0], imgPointList[j][0]),(imgPointList[i][1], imgPointList[j][1]) , 'b-')
    ax1.set_title('Source plane')
    ax2.set_title('Image plane')
    plt.show()





def lm_arctanh(x):
    if x<-1 or x>1:
        print "x should be between -1 and 1"
    return np.log(np.sqrt((1.0+x)/(1.0-x)))



def main():
    return "Nothing to do!"


if __name__=='__main__':
    main()
