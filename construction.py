__author__ = 'cheng109'

import numpy as np
import commons
import scipy.sparse
import scipy.sparse.linalg
## construction for source image vector

def getLensOperator(mappingDict, srcBrightNess):
    dim = len(mappingDict)
    L = np.zeros((dim, 2*dim))
    imgPointList = mappingDict.keys()
    srcPointList = mappingDict.values()
    srcPosition = []
    firstOrderWeightList = []
    secondOrderWeightList = []
    Hsy1 = np.zeros((dim, dim))
    Hsy2 = np.zeros((dim, dim))
    d =  np.zeros(dim)

    normV = [[] for i in range(dim)]

    #srcBrightNess = np.ones(dim)
    varList = []

    RTR = np.zeros((2*dim, 2*dim))

    for i in range(dim):
        imgX, imgY = imgPointList[i]
        srcX, srcY, imgBrightNess, type, var = srcPointList[i]
        srcPosition.append((srcX, srcY))
        varList.append(var)
        d[i] = imgBrightNess

        if type=='o':
            if (imgX, imgY-1) in imgPointList and (imgX, imgY+1) in imgPointList and (imgX+1, imgY) in imgPointList:
                Hsy1[i][i]=1
                Hsy2[i][i]=1
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
                firstOrderWeightList.append((w1Index, w2Index, w3Index, L[i][w1Index], L[i][w2Index],L[i][w3Index]))
            else:
                L[i][i] =1
                firstOrderWeightList.append((i, i, i , 1/3.,1/3.,1/3.))


        if type=='v':
            L[i][i] =1
            firstOrderWeightList.append((i, i, i , 1/3.,1/3.,1/3.))
            if (imgX-1, imgY-1) in imgPointList and (imgX-1, imgY+1) in imgPointList and (imgX+1, imgY-1) in imgPointList and (imgX+1, imgY+1) in imgPointList:
                wAIndex = imgPointList.index((imgX-1, imgY-1))
                wBIndex = imgPointList.index((imgX-1, imgY+1))
                wCIndex = imgPointList.index((imgX+0, imgY+0))
                wDIndex = imgPointList.index((imgX+1, imgY-1))
                wEIndex = imgPointList.index((imgX+1, imgY+1))


                A = (mappingDict[(imgX-1,imgY-1)][0], mappingDict[(imgX-1,imgY-1)][1])
                B = (mappingDict[(imgX-1,imgY+1)][0], mappingDict[(imgX-1,imgY+1)][1])
                C = (mappingDict[(imgX+0,imgY+0)][0], mappingDict[(imgX+0,imgY+0)][1])
                D = (mappingDict[(imgX+1,imgY-1)][0], mappingDict[(imgX+1,imgY-1)][1])
                E = (mappingDict[(imgX+1,imgY+1)][0], mappingDict[(imgX+1,imgY+1)][1])
                Hsy1[i][wAIndex], Hsy1[i][wBIndex], Hsy1[i][wCIndex], Hsy1[i][wDIndex],Hsy1[i][wEIndex],\
                Hsy2[i][wAIndex], Hsy2[i][wBIndex], Hsy2[i][wCIndex], Hsy2[i][wDIndex],Hsy2[i][wEIndex] = commons.getPentWeigth(A, B, C, D, E)


    HsTHs = scipy.sparse.coo_matrix(Hsy1.transpose())*scipy.sparse.coo_matrix(Hsy1) + \
                    scipy.sparse.coo_matrix(Hsy2.transpose())*scipy.sparse.coo_matrix(Hsy2)

    HsTHs = HsTHs.toarray()
    #Hs = np.linalg.cholesky(tt)


    for i in range(dim):
        for j in range(dim):
            RTR[i][j] = HsTHs[i][j]


        # Diagonal covariance matrix:
    C = commons.listToDiagonalMatrix(varList)
    normV = commons.getMeanNorm(normV)

    return srcPosition, srcPointList, commons.sMatrix(L), normV, C, firstOrderWeightList, d, RTR, HsTHs



def getPSFMatrix(psfFileName,filterMatrix,const):
    psfMatrix, _= commons.readFitsImage(psfFileName)
    normalizedPSF = psfMatrix/np.sum(psfMatrix)
    M, N  = const.imgSize
    fullB = np.zeros((M*N, M*N))
    xlim, ylim = normalizedPSF.shape

    #B = np.zeros((const.length, const.length))
    for u in range(xlim):
        for v in range(ylim):
            for h in range(M*N):
                if h+u>=0 and h+u<=M-1 and h+v>=0 and h+v<=N-1:
                    g=(h+u)+(h+v-1)*(M-1)
                    fullB[g][h]=normalizedPSF[u][v]


    B = commons.sMatrix(fullB)*filterMatrix
    return B


def getPenalty(M, r, d, Hs, s, Hphi,phi, lambdaS, lambdaPhi):

    residual = np.outer(M, r)-d
    chiSquare = np.transpose(residual)*np.linalg.inv(C)*residual
    regSource = lambdaS**2*np.linalg.norm(np.outer(Hs, s))**2
    regPot    = lambdaPhi**2*np.linalg.norm(np.outer(Hphi,phi))**2
    penalty = chiSquare+regSource+regPot
    return penalty




def getChiSquare(M, r, d, C, const):

    #M = scipy.sparse.coo_matrix(M)

    r = np.reshape(r, (1,2*const.length))
    d = np.reshape(d, (1,const.length))
    r = scipy.sparse.coo_matrix(np.transpose(r))
    d = scipy.sparse.coo_matrix(np.transpose(d))
    C = scipy.sparse.coo_matrix(C)

    #r = scipy.sparse.linalg.spsolve(M, d)
    res = M*r-d
    resT = res.transpose()
    chiSquare = resT*scipy.sparse.linalg.inv(C)*res
    return chiSquare

def getMatrixDs(normV,indexWeightList, const):
    # Assume the source vector 's' is sorted
    xdim =const.length
    ydim = 2*xdim
    DsMatrix = np.zeros((xdim, ydim))

    for i in range(len(indexWeightList)):

        indexA, indexB, indexC, wA, wB, wC = indexWeightList[i]
        if not normV[i]:
            nA0, nA1, nA2 = normV[indexA]
            nB0, nB1, nB2 = normV[indexB]
            nC0, nC1, nC2 = normV[indexC]
            dSp_y1 = wA*(-nA0/nA2)  + wB*(-nB0/nB2) + wC*(-nC0/nC2)
            dSp_y2 = wA*(-nA1/nA2)  + wB*(-nB1/nB2) + wC*(-nC1/nC2)
        else:
            n0, n1, n2 = normV[i]
            if n2==0:
                dSp_y1=0
                dSp_y2=0
            else:
                dSp_y1 = -n0/n2
                dSp_y2 = -n1/n2
        #print i
        DsMatrix[i][2*i]    = dSp_y1
        DsMatrix[i][2*i+1]  = dSp_y2
    return DsMatrix

# def getMatrixDphi(dPhiGrid, potPosition, D, const):
#     xdim, ydim = const.potSize
#     DphiMatrix = np.zeros((2*xdim, ydim))
#     #dPhiGrid.reshape((xdim, ydim))
#     for i in range(xdim):
#         for j in range(ydim-1):
#             col = i*ydim+j
#             DphiMatrix[2*col][col] = (dPhiGrid[i][j+1]-dPhiGrid[i][j])/const.potRes
#             DphiMatrix[2*col+1][col] = (dPhiGrid[i+1][j]-dPhiGrid[i][j])/const.potRes
#     return DphiMatrix

def getMatrixDphi(const):

    ydim =  const.length

    DphiMatrix = np.zeros((2*ydim, ydim))
    for col in range(ydim-1):
        DphiMatrix[2*col][col] = -1
        DphiMatrix[2*col+1][col] = -1
    return DphiMatrix



def getMatrixH(const):
    ydim = const.length
    HsMatrix = np.zeros((ydim,ydim))



    return HsMatrix

def main():
    return


if __name__ =="__main__":
    main()