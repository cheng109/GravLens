
from astropy.io import fits
from fastell4py import fastell4py
import numpy as np
import matplotlib.pyplot as plt
import pyfits
import pyregion
import scipy.sparse

class Constants():
    def __init__(self, srcSize, imgSize,potSize, srcRes, imgRes, potRes, length):
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
        self.length = length


def sMatrix(m):
    return scipy.sparse.coo_matrix(m)

def filterImage(maskFileName, imageFileName, filterType):
    imgList = []
    imgData, _ = readFitsImage(imageFileName)
    xlim, ylim = imgData.shape
    if filterType=="REG":
        reg = open(maskFileName, 'r').read()

        hdulist = pyfits.open(imageFileName)
        maskData = pyregion.parse(reg).get_mask(hdu=hdulist[0])
        maskData = writeMaskFile(maskData, 'mask.fits')
    if filterType=="FITS":
        maskData, _= readFitsImage(maskFileName)
        assert (maskData.shape==imgData.shape), "Mask shape does not match image shape!"

    if filterType=="NONE":
        maskData = np.ones((xlim, ylim))

    for i in range(xlim):
        for j in range(ylim):
            if maskData[i][j]!=0:
                if (i+j)%2 ==1:
                    type = 'v'
                else:
                    type = 'o'
                imgList.append((i,j, imgData[i][j], type))

    return imgList

def getFilterMatrix(imgList, const):
    xlim,ylim = const.imgSize
    filter_1 = np.zeros((xlim*ylim,len(imgList)))
    filter_2 = np.zeros((len(imgList),xlim*ylim))

    for k in range(len(imgList)):
        i, j, _, _ = imgList[k]
        filter_1[i*ylim+j][k] = 1

    for k in range(len(imgList)):
        i, j, _, _ = imgList[k]
        filter_2[k][i*ylim+j] = 1

    return filter_1, filter_2

def writeMaskFile(maskData, outputName):
    xlim, ylim = maskData.shape
    mask = np.zeros((xlim, ylim))
    for i in range(xlim):
        for j in range(ylim):
            if maskData[i][j]==True:
                mask[i][j] =1
            #else:
            #[i][j] = 0
    writeFitsImage(mask, outputName)
    return mask

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
        if x>0 and x<const.srcSize[0]-1 and y>0 and  y<const.srcSize[1]-1:
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

def getPentWeigth(A, B, C, D, E):
    # C is the center point
    Q = getLinerInterpolate(A, B,C, direction='x')
    P = getLinerInterpolate(D, E, C, direction='x')

    M = getLinerInterpolate(B, E, C, direction='y')
    N = getLinerInterpolate(A, D, C, direction='y')

    XweightA = dist(Q, B)/(dist(C,Q)*dist(A,B))
    XweightB = dist(Q, A)/(dist(C,Q)*dist(A,B))
    XweightC = -(1/dist(C,P)+1/dist(C,Q))
    XweightD = dist(P, E)/(dist(C, P)*dist(D, E))
    XweightE = dist(P, D)/(dist(C, P)*dist(D, E))

    YweightA = dist(N, D)/(dist(C, N)*dist(A,D))
    YweightB = dist(M, E)/(dist(C,M)*dist(B,E))
    YweightC = -(1/dist(C,N)+1/dist(C,M))
    YweightD = dist(A, N)/(dist(C, N)*dist(A,D))
    YweightE = dist(B, M)/(dist(C, M)*dist(B,E))

    return XweightA,XweightB,XweightC,XweightD,XweightE, YweightA,YweightB,YweightC,YweightD,YweightE

def dist(A, B):
    return np.sqrt((A[0]-B[0])**2+(A[1]-B[1])**2)

def getLinerInterpolate(A, B, C, direction):
    xA,yA = A
    xB,yB = B
    xC,yC = C

    if abs(xA-xB) < 10e-8:
        P = ((xA+xB)*0.5, yC)
    else:
        a = float(yA-yB)/(xA-xB)
        b = yA - a*xA
        if direction=='x':
            P = ((yC-b)/a, yC)
        if direction=='y':
            P = (xC, a*xC+b)
    return P


def getMeanNorm(normV):
    meanNorm = []
    for i in range(len(normV)):
        if len(normV[i])==0:
            meanNorm.append((0,0,0))
        else:
            sumN0,sumN1, sumN2, counter = 0, 0, 0, 0
            for j in range(len(normV[i])):
                n0, n1, n2 = normV[i][j]
                sumN0 += n0
                sumN1 += n1
                sumN2 += n2
                counter += 1
            meanNorm.append((sumN0/counter, sumN1/counter, sumN2/counter))
    return meanNorm

def getNormVectors(p1,p2, p3):
    p1x, p1y, p1z = p1
    p2x, p2y, p2z = p2
    p3x, p3y, p3z = p3

    nx = (p3y-p2y)*(p2z-p1z)-(p3z-p2z)*(p2y-p1y)
    ny = (p3z-p2z)*(p2x-p1x)-(p3x-p2x)*(p2z-p1z)
    nz = (p3x-p2x)*(p2y-p1y)-(p3y-p2y)*(p2x-p1x)
    return nx, ny, nz


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
        srcX, srcY, _ , type, _ = srcPointList[i]

        if type=='v' and (imgX, imgY+2) in mappingDict and (imgX+1, imgY+1) in mappingDict:

            ax1.plot((srcX, mappingDict[(imgX, imgY+2)][0]), (srcY, mappingDict[(imgX, imgY+2)][1]), 'b-')
            ax1.plot((srcX, mappingDict[(imgX+1, imgY+1)][0]), (srcY, mappingDict[(imgX+1, imgY+1)][1]), 'b-')
            ax1.plot((mappingDict[(imgX, imgY+2)][0], mappingDict[(imgX+1, imgY+1)][0]), (mappingDict[(imgX, imgY+2)][1], mappingDict[(imgX+1, imgY+1)][1]), 'b-')
        # plot  the lensed image plane grid
        for j in np.arange(i+1, len(imgPointList), 1):
            if (imgPointList[i][0]==imgPointList[j][0] and abs(imgPointList[i][1]-imgPointList[j][1])==1) or (imgPointList[i][1]==imgPointList[j][1] and abs(imgPointList[i][0]-imgPointList[j][0])==1):
                #ax1.plot((srcPointList[i][0], srcPointList[j][0]),(srcPointList[i][1], srcPointList[j][1]) , 'b-')
                ax2.plot((imgPointList[i][0], imgPointList[j][0]),(imgPointList[i][1], imgPointList[j][1]) , 'b-')
        #plot the 'vertex' and 'ohters'
        if type=='v':
            ax2.plot(imgX, imgY, 'ro')
            ax1.plot(srcX, srcY, 'ro')
        else:
            ax1.plot(srcX, srcY, 'wo')
            ax2.plot(imgX, imgY, 'wo')

    ax1.set_title('Source plane')
    ax2.set_title('Image plane')

    ax1.set_xlim([0, const.srcSize[0]])
    ax1.set_ylim([0, const.srcSize[1]])

    plt.show()



def disp(data):
    plt.imshow(data, origin="lower", interpolation="nearest")
    plt.show()

def lm_arctanh(x):
    if x<-1 or x>1:
        print "x should be between -1 and 1"
    return np.log(np.sqrt((1.0+x)/(1.0-x)))

def main():

    #filter= createGirdFilter(50, 50 )
    #plt.imshow(filter, interpolation="nearest")
    #plt.show()
    A = (-1, -2)
    B = (-1, 1)
    C = (0, 0)
    D = (1, -1)
    E = (1, 1)


    print getPentWeigth(A, B, C, D, E)

    getLinerInterpolate(A,D,C,direction='y')
    return "Nothing to do!"


if __name__=='__main__':
    main()
