import commons
from fastell4py import fastell4py
import random
import numpy as np
import deflection
import construction
import matplotlib.pyplot as plt


def inverseMapping(srcData,imgData, mappingImage, lens, const):
    #### mappingDict={'imageGrid': srcGrid,  'imageGrid':srcGrid, .....}
    mappingDict = {}
    maskFileName = 'test_images/mask.fits'
    maskData, mVector = commons.readFitsImage(maskFileName)




    for i in range(mappingImage.shape[1]):
        for j in range(mappingImage.shape[0]):
            if maskData[i][j]==1:
                _, _, srcX, srcY  = deflection.getDeflection(i, j, 'SIE',const)
                if (i+j)%2 ==1:
                    type = 'v'
                else:
                    type = 'o'
                mappingDict[(i,j)]= (srcX, srcY, imgData[i][j], type)
                if srcX<const.srcSize[0] and srcY<const.srcSize[1] and srcX >0 and srcY>0:
                    mappingImage[i][j]= srcData[srcX][srcY]
    #commons.plotMappingDict(mappingDict, const)

    return mappingDict, mappingImage




def main():
    dirName= "test_images/"
    imgData, dVector = commons.readFitsImage(dirName + "model_img.fits")

    varData, vVector = commons.readFitsImage(dirName+ "var.fits.gz")
    B, _= commons.readFitsImage(dirName+ "psf.fits.gz")
    srcData = np.zeros((80,80))
    const = commons.Constants(srcData.shape,imgData.shape, potSize=potSize,srcRes=0.02, imgRes=0.05, potRes=0.05)

    #plt.imshow(np.reshape(dVector, (imgData.shape[0], imgData.shape[1])),origin="lower")
    #plt.imshow(B, origin="lower", interpolation="nearest")

    mappingDict, mappingImage = inverseMapping(srcData, imgData,  mappingImage=imgData, lens=0, const=const)
    srcPosition, srcBrightNess, sXDev, sYDev = construction.getLensOperator(mappingDict)
    srcMap = commons.recoverSource(srcPosition, srcBrightNess, const)
    plt.imshow(srcMap, origin="lower", interpolation="nearest")
    plt.show()
    #commons.writeFitsImage(srcMap, "model_test.fits")


if __name__=='__main__':
    main()
