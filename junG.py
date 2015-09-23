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


def getMappingDict(imgList,varList, model, const):
    mappingDict = {}

    for index in range(len(imgList)):
        i, j, imgBrightNess, type = imgList[index]
        _, _, srcX, srcY = deflection.getDeflection(i, j, model=model, const=const)
        mappingDict[(i,j)]= (srcX, srcY, imgBrightNess, type, varList[index][2])
    commons.plotMappingDict(mappingDict, const)

    return mappingDict


def prepare():

    maskFileName = 'test_images/mask.fits'
    imgFileName = 'test_images/image.fits.gz'
    varFileName = 'test_images/var.fits.gz'
    psfFileName = 'test_images/psf.fits.gz'
    const = commons.Constants(srcSize=(80,80),imgSize=(122,122), potSize=(122,122),srcRes=0.02, imgRes=0.05, potRes=0.05)
    imgList = commons.filterImage(maskFileName=maskFileName,imageFileName=imgFileName, filterType="FITS")
    varList = commons.filterImage(maskFileName=maskFileName,imageFileName=varFileName, filterType="FITS")

    s = np.zeros(const.imgSize[0]*const.imgSize[1])
    B = construction.getPSFMatrix(psfFileName, const)

    return s, B, const,  imgList, varList


def procedure():
    # initialization:

    s, B, const,  imgList, varList = prepare()
    mappingDict = getMappingDict(imgList,varList, model='SIE', const=const)
    srcPosition, srcPointList, L, normV, C, indexWeightList, d = construction.getLensOperator(mappingDict)
    print normV
    Ds = construction.getMatrixDs(normV,indexWeightList, const)

    M =1
    R =1
    temp = np.transpose(M)*np.linalg.inv(C)
    l = len(d)
    A = np.zeros((len(d),len(d)))
    for i in range(l):
        A[i][i] =1

    s = linalg.solve(scipy.matrix(A), d)
    print s.shape
    #s = linalg.spsolve(lil_matrix(np.outer(temp,M) + R).tocsr(), np.outer(temp,d))
    print "s= ", s

    srcMap = commons.pixelizeSource(srcPosition, s, const)
    print srcMap
    commons.writeFitsImage(srcMap, "model_test.fits")
    plt.imshow(srcMap, origin="lower", interpolation="nearest")
    plt.show()
    return


def main():
    procedure()

if __name__=='__main__':
    main()
