import commons
from fastell4py import fastell4py
import random
import scipy
import numpy as np
import deflection
import construction
import matplotlib.pyplot as plt
#from scipy import linalg
#from scipy.sparse import linalg
from scipy import linalg
from scipy.sparse import lil_matrix
from numpy.linalg import norm
import models

def getMappingDict(imgList,varList, MODEL , const):
    mappingDict = {}

    for index in range(len(imgList)):
        i, j, imgBrightNess, type = imgList[index]
        _, _, srcX, srcY = deflection.getDeflection(i+0.5, j+0.5, MODEL, const=const)
        mappingDict[(i,j)]= (srcX, srcY, imgBrightNess, type, varList[index][2])
    #commons.plotMappingDict(mappingDict, const)

    return mappingDict


def prepare():

    dir = 'example/'

    maskFileName = dir + 'mask_2.reg'
    imgFileName =  dir + 'jun_image.fits'
    varFileName =  dir + 'jun_var.fits'
    psfFileName =  dir + 'jun_psf.fits'

    maskFilterType = 'REG'

    imgList = commons.filterImage(maskFileName=maskFileName,imageFileName=imgFileName, filterType=maskFilterType)
    varList = commons.filterImage(maskFileName=maskFileName,imageFileName=varFileName, filterType=maskFilterType)



    imgSize = commons.getImageSize(imgFileName)

    const = commons.Constants(srcSize=imgSize,imgSize=imgSize, potSize=imgSize,srcRes=0.3, imgRes=0.3, potRes=0.3, length=len(imgList))
    s = np.zeros(len(imgList))
    phi = np.zeros(len(imgList))
    B = construction.getPSFMatrix(psfFileName, imgList, const)
    lambdaS = 3


    return s, B, const,  imgList, varList, phi, lambdaS


def procedure():
    # initialization:

    s, B, const,  imgList, varList, phi, lambdaS = prepare()

    LM_SIE = models.PTMASS(critRad=6.16, centerX=0.0, centerY=5.0)
    #LM_SIE = models.SIE(critRad=6.16, axisRatio=0.7, pa=0.0, centerX=0.0, centerY=0.0)

    mappingDict = getMappingDict(imgList,varList, LM_SIE, const=const)

    imgPointList = mappingDict.keys()
    srcPointList = mappingDict.values()
    dim = len(imgPointList)
    srcPosition = []
    s = []
    r = np.zeros((1, 2*dim))
    for i in range(dim):
        x,y, bri, _, _ = srcPointList[i]

        srcPosition.append((x,y))
        s.append(bri)
        r[0][i] = bri
    srcMap = commons.pixelizeSource(srcPosition, s, const)


    commons.writeFitsImage(srcMap, "example/outputSrc.fits")
    plt.imshow(srcMap, origin="lower", interpolation="nearest")
    plt.colorbar()
    plt.show()


    #return
    chi2List = []
    #critRadRange = np.arange(5.9, 6.3, 0.02)
    #qRange = np.arange(0.1, 0.9, 0.05)
    #paRange = np.arange(0, 20, 10)
    centerYrange= np.arange(5, 15, 1)
    index = []
    i = 0
    # = []
    #for critRad in critRadRange:
    # r = np.reshape(r, (1,2*const.length))
    # r = scipy.sparse.coo_matrix(np.transpose(r))

    for yshift in centerYrange:

        i = i+1
        #LM_SIE.critRad = critRad
        #LM_SIE.axisRatio = q
        #LM_SIE.pa = pa
        LM_SIE.centerY = yshift
        mappingDict = getMappingDict(imgList,varList, MODEL = LM_SIE, const=const)


        srcPosition, srcPointList, L, normV, C, indexWeightList, d, RTR , HsTHs, imgPointList= construction.getLensOperator(mappingDict, s)

        chi2, res = construction.getChiSquare(M = L, r = r, d =d,  C=C, const=const)
        chi2 = chi2.toarray()[0][0]
        res = res.toarray()

        # image = construction.buildModelImage(imgPointList, res, const)
        # plt.imshow(image, origin="lower", interpolation="nearest")
        # plt.colorbar()
        # plt.show()
        print  yshift,  chi2


        index.append(yshift)
        chi2List.append(chi2)





    plt.plot(index, chi2List, '*-')
    plt.xlabel("Parameter")
    plt.ylabel("Chi2")
    plt.show()

    #critRadRange = np.arange(1.3, 4.0, 0.1)

    # for i in range(1):
    #     mappingDict = getMappingDict(imgList,varList, model='SIE', critRad = 3.0, const=const)
    #     srcPosition, srcPointList, L, normV, C, indexWeightList, d, RTR, HsTHs = construction.getLensOperator(mappingDict, s)
    #     Ds = construction.getMatrixDs(normV,indexWeightList, const)
    #     Dphi =construction.getMatrixDphi(const)
    #
    #
    #     D = scipy.sparse.coo_matrix(Ds)*scipy.sparse.coo_matrix(Dphi)
    #
    #
    #     M =scipy.sparse.hstack([L,D])
    #
    #     r = np.zeros(2*const.length)
    #     for i in range(const.length):
    #         r[i] = s[i]
    #
    #
    #     r =scipy.sparse.vstack(r)
    #     M_t = scipy.sparse.coo_matrix.transpose(M)
    #
    #     C = scipy.sparse.coo_matrix(C)
    #     RTR = scipy.sparse.coo_matrix(RTR)
    #
    #     C_inv = scipy.sparse.linalg.inv(C)
    #     b = M_t*C_inv*d
    #
    #     A = M_t*C_inv*M+RTR
    #     r = scipy.sparse.linalg.spsolve(A, b)
    #
    #     s = np.zeros((const.length,1))
    #     for i in range(const.length):
    #         s[i][0] = r[i]
    #     s_t = scipy.sparse.coo_matrix.transpose(scipy.sparse.coo_matrix(s))
    #     chi2 = construction.getChiSquare(M = M, r = r, d =d,  C=C)
    #     G = chi2 + s_t*scipy.sparse.coo_matrix(HsTHs)*s
    #     print i, chi2, G
    # srcMap = commons.pixelizeSource(srcPosition, s, const)
    # print srcMap
    # commons.writeFitsImage(srcMap, "model_test.fits")
    # plt.imshow(srcMap, origin="lower", interpolation="nearest")
    # plt.show()



    return


def main():
    procedure()

if __name__=='__main__':
    main()
