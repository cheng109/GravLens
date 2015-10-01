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


def getMappingDict(imgList,varList, model, critRad , const):
    mappingDict = {}

    for index in range(len(imgList)):
        i, j, imgBrightNess, type = imgList[index]
        _, _, srcX, srcY = deflection.getDeflection(i, j, model=model, critRad = critRad,const=const)
        mappingDict[(i,j)]= (srcX, srcY, imgBrightNess, type, varList[index][2])
    #commons.plotMappingDict(mappingDict, const)

    return mappingDict


def prepare():

    maskFileName = 'test_images/mask.fits'
    imgFileName = 'test_images/image.fits.gz'
    varFileName = 'test_images/var.fits.gz'
    psfFileName = 'test_images/psf.fits.gz'


    imgList = commons.filterImage(maskFileName=maskFileName,imageFileName=imgFileName, filterType="FITS")
    varList = commons.filterImage(maskFileName=maskFileName,imageFileName=varFileName, filterType="FITS")

    const = commons.Constants(srcSize=(80,80),imgSize=(122,122), potSize=(122,122),srcRes=0.02, imgRes=0.05, potRes=0.05, length=len(imgList))

    s = np.zeros(len(imgList))
    phi = np.zeros(len(imgList))
    B = construction.getPSFMatrix(psfFileName, const)

    return s, B, const,  imgList, varList, phi


def procedure():
    # initialization:

    s, B, const,  imgList, varList, phi = prepare()


    chi2List = []
    # critRadRange = np.arange(0.5, 2.0, 0.02)
    # index = []
    # i = 0
    # for critRad in critRadRange:
    #     i = i+1
    #     mappingDict = getMappingDict(imgList,varList, model='SIE', critRad = critRad, const=const)
    #     srcPosition, srcPointList, L, normV, C, indexWeightList, d, RTR = construction.getLensOperator(mappingDict)
    #
    #
    #     chi2 = construction.getChiSquare(M = L, r = d, d =d,  C=C)
    #     print  critRad, chi2.shape, chi2[0][0]
    #     index.append(critRad)
    #     chi2List.append(chi2[0][0])
    #
    # plt.plot(index, chi2List, '*-')
    # plt.show()

    critRadRange = np.arange(1.3, 4.0, 0.1)

    for i in range(1):
        mappingDict = getMappingDict(imgList,varList, model='SIE', critRad = 3.0, const=const)
        srcPosition, srcPointList, L, normV, C, indexWeightList, d, RTR, HsTHs = construction.getLensOperator(mappingDict, s)
        Ds = construction.getMatrixDs(normV,indexWeightList, const)
        Dphi =construction.getMatrixDphi(const)


        D = scipy.sparse.coo_matrix(Ds)*scipy.sparse.coo_matrix(Dphi)


        M =scipy.sparse.hstack([L,D])

        r = np.zeros(2*const.length)
        for i in range(const.length):
            r[i] = s[i]


        r =scipy.sparse.vstack(r)
        M_t = scipy.sparse.coo_matrix.transpose(M)

        C = scipy.sparse.coo_matrix(C)
        RTR = scipy.sparse.coo_matrix(RTR)

        C_inv = scipy.sparse.linalg.inv(C)
        b = M_t*C_inv*d

        A = M_t*C_inv*M+RTR
        r = scipy.sparse.linalg.spsolve(A, b)

        s = np.zeros((const.length,1))
        for i in range(const.length):
            s[i][0] = r[i]
        s_t = scipy.sparse.coo_matrix.transpose(scipy.sparse.coo_matrix(s))
        chi2 = construction.getChiSquare(M = M, r = r, d =d,  C=C)
        G = chi2 + s_t*scipy.sparse.coo_matrix(HsTHs)*s
        print i, chi2, G
    # srcMap = commons.pixelizeSource(srcPosition, s, const)
    # print srcMap
    # commons.writeFitsImage(srcMap, "model_test.fits")
    # plt.imshow(srcMap, origin="lower", interpolation="nearest")
    # plt.show()

    s_stable =s
    critList = []
    GList = []
    for critRad in critRadRange:
        s = s_stable
        mappingDict = getMappingDict(imgList,varList, model='SIE', critRad = critRad, const=const)
        srcPosition, srcPointList, L, normV, C, indexWeightList, d, RTR, HsTHs = construction.getLensOperator(mappingDict, s)
        Ds = construction.getMatrixDs(normV,indexWeightList, const)
        Dphi =construction.getMatrixDphi(const)


        D = scipy.sparse.coo_matrix(Ds)*scipy.sparse.coo_matrix(Dphi)


        M =scipy.sparse.hstack([L,D])

        r = np.zeros(2*const.length)
        for i in range(const.length):
            r[i] = s[i]


        r =scipy.sparse.vstack(r)
        M_t = scipy.sparse.coo_matrix.transpose(M)

        C = scipy.sparse.coo_matrix(C)
        RTR = scipy.sparse.coo_matrix(RTR)

        C_inv = scipy.sparse.linalg.inv(C)
        b = M_t*C_inv*d

        A = M_t*C_inv*M+RTR
        r = scipy.sparse.linalg.spsolve(A, b)

        s = np.zeros((const.length,1))
        for i in range(const.length):
            s[i][0] = r[i]
        s_t = scipy.sparse.coo_matrix.transpose(scipy.sparse.coo_matrix(s))
        chi2 = construction.getChiSquare(M = M, r = r, d =d,  C=C)
        G = chi2 + s_t*scipy.sparse.coo_matrix(HsTHs)*s
        print critRad, chi2[0][0], G[0][0]
        critList.append(critRad)
        GList.append(G[0][0])
    plt.plot(critList, GList, '*-')

    return


def main():
    procedure()

if __name__=='__main__':
    main()
