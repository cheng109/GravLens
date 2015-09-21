
from astropy.io import fits
from fastell4py import fastell4py
import numpy as np
import matplotlib.pyplot as plt


class Constants():
    def __init__(self, srcSize, imgSize,potSize, srcRes, imgRes, potRes):
        self.srcSize = srcSize
        self.imgSize = imgSize
        self.potSize = potSize
        self.srcRes = srcRes
        self.imgRes = imgRes
        self.potRes = potRes


        self.srcXCenter = srcSize[0]/2.0
        self.srcYCenter = srcSize[1]/2.0
        self.imgXCenter = imgSize[0]/2.0
        self.imgYCenter = imgSize[1]/2.0



def readFitsImage(imageName):
    hdulist = fits.open(imageName)
    data =  hdulist[0].data
    hdulist.close()
    vector = np.reshape(data, data.shape[0]*data.shape[1])
    return data, vector

def writeFitsImage(data, outputName):
    hdu = fits.PrimaryHDU(data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(outputName, clobber=True)

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

    maskData, mVector = readFitsImage(maskFileName)
    newMappingDict = {}
    for i in range(maskData.shape[0]):
        for j in range(maskData.shape[1]):
            if maskData[i][j]>0.5:
                newMappingDict[(i,j)]=mappingDict[(i,j)]
    return newMappingDict


def createGirdFilter(xlen, ylen):

    filter = np.zeros([xlen,ylen])
    for i in range(xlen):
        for j in range(ylen):
            if (i+j)%2==1:
                filter[i][j]=1
    return filter

def pixelizeSource(srcPosition, srcBrightNess , const):

    srcMap = np.zeros(const.srcSize)
    for i in range(len(srcPosition)):
        x, y = srcPosition[i]
        if x<const.srcSize[0]-1 and y<const.srcSize[1]-1:
            srcMap[int(x)][int(y)] +=  srcBrightNess[i]

    # return a pixelized source map.
    return srcMap

def getTriWeight(A,B,C, P):
    def area(a, b, c):
        def distance(p1, p2):
            return np.hypot(p1[0]-p2[0], p1[1]-p2[1])

        side_a = distance(a, b)
        side_b = distance(b, c)
        side_c = distance(c, a)
        s = 0.5 * ( side_a + side_b + side_c)
        return np.sqrt(s * (s - side_a) * (s - side_b) * (s - side_c))
    areaA = area(P, B, C)
    areaB = area(P, A, C)
    areaC = area(P, A, B)
    S = areaA+areaB + areaC
    return areaA/S, areaB/S, areaC/S



def getMeanNorm(normV):
    meanNorm = []
    for i in range(len(normV)):
        if len(normV[i])==0:
            meanNorm.append((0,0,0))
        sumN0,sumN1, sumN2, counter = 0, 0, 0, 0
        for (n0, n1, n2) in normV[i]:
            sumN0 += n0
            sumN1 += n1
            sumN2 += n2
            counter += 1
        meanNorm.append((sumN0/counter, sumN1/counter, sumN2/counter))
    return meanNorm


### Build lens operator



def getNormVectors(A, B, O):

    a1 = A[0] - O[0]
    a2 = A[1] - O[1]
    a3 = A[2] - O[2]
    b1 = B[0] - O[0]
    b2 = B[1] - O[1]
    b3 = B[2] - O[2]

    n0 = a2*b3-a3*b2
    n1 = -(a1*b3-a3*b1)
    n2 = a1*b2-a2*b1
    return n0, n1, n2


def listToDiagonalMatrix(l):
    dim = len(l)
    D = np.zeros((dim,dim))
    for i in range(dim):
        D[i][i]= l[i]
    return D


def plotMappingDict(mappingDict,const):
   #### mappingDict={'imageGrid': srcGrid,  'imageGrid':srcGrid, .....}
    f, (ax1, ax2) = plt.subplots(1, 2) #, sharex=True, sharey=True )

    imgPointList = mappingDict.keys()
    srcPointList = mappingDict.values()




    for i in range(len(imgPointList)):
        imgX, imgY = imgPointList[i]
        srcX, srcY, _ , type = srcPointList[i]

        if type=='v' and (imgX, imgY+2) in mappingDict and (imgX+1, imgY+1) in mappingDict:

            ax1.plot((srcX, mappingDict[(imgX, imgY+2)][0]), (srcY, mappingDict[(imgX, imgY+2)][1]), 'b-')
            ax1.plot((srcX, mappingDict[(imgX+1, imgY+1)][0]), (srcY, mappingDict[(imgX+1, imgY+1)][1]), 'b-')
            ax1.plot((mappingDict[(imgX, imgY+2)][0], mappingDict[(imgX+1, imgY+1)][0]), (mappingDict[(imgX, imgY+2)][1], mappingDict[(imgX+1, imgY+1)][1]), 'b-')
        # plot  the lensed image plane grid
        for j in np.arange(i+1, len(imgPointList), 1):
            if (imgPointList[i][0]==imgPointList[j][0] and abs(imgPointList[i][1]-imgPointList[j][1])==1) or (imgPointList[i][1]==imgPointList[j][1] and abs(imgPointList[i][0]-imgPointList[j][0])==1):
                #ax1.plot((srcPointList[i][0], srcPointList[j][0]),(srcPointList[i][1], srcPointList[j][1]) , 'b-')
                ax2.plot((imgPointList[i][0], imgPointList[j][0]),(imgPointList[i][1], imgPointList[j][1]) , 'b-')
        # plot the 'vertex' and 'ohters'
        if type=='v':
            ax2.plot(imgX, imgY, 'ro')
            ax1.plot(srcX, srcY, 'ro')
        else:
            ax1.plot(srcX, srcY, 'wo')
            ax2.plot(imgX, imgY, 'wo')

    ax1.set_title('Source plane')
    ax2.set_title('Image plane')

    plt.show()





def lm_arctanh(x):
    if x<-1 or x>1:
        print "x should be between -1 and 1"
    return np.log(np.sqrt((1.0+x)/(1.0-x)))

def main():

    filter= createGirdFilter(50, 50 )
    plt.imshow(filter, interpolation="nearest")
    plt.show()
    return "Nothing to do!"


if __name__=='__main__':
    main()
