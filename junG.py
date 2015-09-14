import commons
from fastell4py import fastell4py
import random
import numpy as np
import deflection

imgRes = 0.05
srcRes = 0.15
srcDimX=50
srcDimY=50

def getDeflectionAngle(pixImageX, pixImageY):

    return pixSrcX, pixSrcY

def inverseMapping(srcData,imgData, mappingImage, lens, const):
    srcSize = srcData.shape
    imgSize = imgData.shape
    mappingDict = {}
    # i, j are indices in image plane
    # I, J are indices in source plane
    for i in range(mappingImage.shape[1]):
        for j in range(mappingImage.shape[0]):
            _, _, I, J = deflection.getDeflection(i, j, 'SPEMD',const)
            mappingDict[(i,j)]= (I, J)
            if I<const.srcSize[0] and J<const.srcSize[1] and I >0 and J>0:
                mappingImage[i][j]= srcData[I][J]

    #### mappingDict={'imageGrid': srcGrid,  'imageGrid':srcGrid, .....}
    maskFileName = 'test_images/mask.fits'
    newMappingDict = commons.applyMask(maskFileName, mappingDict)
    commons.plotMappingDict(newMappingDict, const)



    return newMappingDict, mappingImage




def main():
    dirName= "test_images/"
    imgData = commons.readFitsImage(dirName + "model_img.fits")
    srcData = commons.readFitsImage(dirName + "model_src.fits")
    varData = commons.readFitsImage(dirName+ "var.fits.gz")
    psfData = commons.readFitsImage(dirName+ "psf.fits.gz")
    const = commons.Constants(srcData.shape,imgData.shape,0.15, 0.05)

    #getModelImage()
    mappingDict, mappingImage = inverseMapping(srcData, imgData,  mappingImage=imgData, lens=0, const=const)
    commons.writeFitsImage(mappingImage, "model_test.fits")





if __name__=='__main__':
    main()
