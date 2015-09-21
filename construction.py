__author__ = 'cheng109'

import numpy as np
import commons
## construction for source image vector



def getLensOperator(mappingDict):
    dim = len(mappingDict)
    L = np.zeros((dim, dim))
    imgPointList = mappingDict.keys()
    srcPointList = mappingDict.values()
    srcPosition = []
    d =  np.zeros(dim)

    normV = [[] for i in range(dim)]

    srcBrightNess = np.ones(dim)
    varList = []
    sXDev = np.zeros(dim)
    sYDev = np.zeros(dim)
    for i in range(dim):
        imgX, imgY = imgPointList[i]
        srcX, srcY, imgBrightNess, type, var = srcPointList[i]
        srcPosition.append((srcX, srcY))
        varList.append(var)
        d[i] = imgBrightNess
        if (imgX, imgY+1) in imgPointList:
            x1, y1 = srcX, srcY
            x2, y2 = (mappingDict[(imgX,imgY+1)][0], mappingDict[(imgX,imgY+1)][1])
            index = imgPointList.index((imgX, imgY+1))
            #sXDev[i] = (srcBrightNess[index]-srcBrightNess[i])/(x2-x1)
            #sYDev[i] = (srcBrightNess[index]-srcBrightNess[i])/(y2-y1)

        if type=='o' and (imgX, imgY-1) in imgPointList and (imgX, imgY+1) in imgPointList and (imgX+1, imgY) in imgPointList:
            w1Index = imgPointList.index((imgX, imgY-1))
            w2Index = imgPointList.index((imgX, imgY+1))
            w3Index = imgPointList.index((imgX+1, imgY))
            A = (mappingDict[(imgX,imgY-1)][0], mappingDict[(imgX,imgY-1)][1])
            B = (mappingDict[(imgX,imgY+1)][0], mappingDict[(imgX,imgY+1)][1])
            C = (mappingDict[(imgX+1,imgY)][0], mappingDict[(imgX+1,imgY)][1])

            L[i][w1Index], L[i][w2Index],L[i][w3Index] = commons.getTriWeight(A,B,C, (srcX, srcY))
            # update the normV
            An = (A[0], A[1], srcBrightNess[w1Index])
            Bn = (B[0], B[1], srcBrightNess[w2Index])
            Cn = (C[0], C[1], srcBrightNess[w3Index])

            n0, n1, n2 = commons.getNormVectors(An, Bn, Cn)
            normV[w1Index].append((n0, n1, n2))
            normV[w2Index].append((n0, n1, n2))
            normV[w3Index].append((n0, n1, n2))
        else:
            L[i][i] =1

        # Diagonal covariance matrix:
        C = commons.listToDiagonalMatrix(varList)



    return srcPosition, srcPointList, L, normV, C



def getPSFMatrix(psfFileName):
    B, _= commons.readFitsImage(psfFileName)
    return B

def getVarianceMatrix(varFileName):


    return


def getMatrixDs():
    # Assume the source vector 's' is sorted


    xdim, ydim = const.potSize
    DsMatrix = np.zeros((xdim, ydim))
    for i in range(xdim):
        for j in range(ydim):

            DsMatrix[i][2*i] = sXDev[i]
            DsMatrix[i][2*i+1]=sYDev[i]


    return DsMatrix


def initialDphi(d, imgPosition):
    dPhi = d
    potPosition = imgPosition

def getMatrixDphi(dPhiGrid, potPosition, D, const):
    xdim, ydim = const.potSize
    DphiMatrix = np.zeros((2*xdim, ydim))
    #dPhiGrid.reshape((xdim, ydim))
    for i in range(xdim):
        for j in range(ydim-1):
            col = i*ydim+j
            DphiMatrix[2*col][col] = (dPhiGrid[i][j+1]-dPhiGrid[i][j])/const.potRes
            DphiMatrix[2*col+1][col] = (dPhiGrid[i+1][j]-dPhiGrid[i][j])/const.potRes
    return DphiMatrix

def main():

    return


if __name__ =="__main__":
    main()