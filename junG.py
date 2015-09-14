import commons
from fastell4py import fastell4py
import random
import numpy as np
import deflection

def inverseMapping(srcData,imgData, mappingImage, lens, const):
    srcSize = srcData.shape
    imgSize = imgData.shape
    mappingDict = {}

    for i in range(mappingImage.shape[1]):
        for j in range(mappingImage.shape[0]):
            _, _, srcX, srcY  = deflection.getDeflection(i, j, 'SIE',const)
            mappingDict[(i,j)]= (srcX, srcY, imgData[i][j])
            if srcX<const.srcSize[0] and srcY<const.srcSize[1] and srcX >0 and srcY>0:
                mappingImage[i][j]= srcData[srcX][srcY]

    #### mappingDict={'imageGrid': srcGrid,  'imageGrid':srcGrid, .....}
    maskFileName = 'test_images/mask.fits'
    newMappingDict = commons.applyMask(maskFileName, mappingDict)
    #commons.plotMappingDict(newMappingDict, const)



    return newMappingDict, mappingImage


def recoverSource( newMappingDict, const):
    srcBrightness = np.zeros(const.srcSize)

    for key, value in newMappingDict.iteritems():
        if value[0]> 0 and value[1] >0 :
            srcBrightness[int(value[0])][int(value[1])] += value[2]
    return srcBrightness





def main():
    dirName= "test_images/"
    imgData = commons.readFitsImage(dirName + "model_img.fits")
    srcData = commons.readFitsImage(dirName + "model_src.fits")
    varData = commons.readFitsImage(dirName+ "var.fits.gz")
    psfData = commons.readFitsImage(dirName+ "psf.fits.gz")
    const = commons.Constants(srcData.shape,imgData.shape,srcRes=0.05, imgRes=0.05)

    #getModelImage()
    mappingDict, mappingImage = inverseMapping(srcData, imgData,  mappingImage=imgData, lens=0, const=const)
    srcBrightness = recoverSource(mappingDict, const)
    commons.writeFitsImage(srcBrightness, "model_test.fits")





if __name__=='__main__':
    main()
