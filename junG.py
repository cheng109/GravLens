import commons
from fastell4py import fastell4py
import random
import numpy as np
import matplotlib.pyplot as plt

imgRes = 0.05
srcRes = 0.15
srcDimX=50
srcDimY=50
def getModelImage():
    pixelratio = 3
    dirName= "test_images/"
    srcData = commons.readFitsImage(dirName + "model_src.fits")
    dimX = srcData.shape[0]
    dimY = srcData.shape[1]

    model = 0
    return model


def getDeflectionAngle(pixImageX, pixImageY):

    critRad = 1.20347*10
    ellipticity = 0.177825
    arat = 1-ellipticity
    PA = 67.5419
    gam = 0.58
    coreRad = 0.05
    q = 0.5*critRad*pow((2-2*gam)/arat, gam)
    s = coreRad**2

    angleImageX = (pixImageX-61)*imgRes
    angleImageY = (pixImageY-61)*imgRes

    deflx, defly = fastell4py.fastelldefl(angleImageX, angleImageY, q, gam, arat, s)
    pixSrcX=int( (angleImageX - deflx)/srcRes + srcDimX/2)
    pixSrcY= int ((angleImageY - defly)/srcRes + srcDimY/2)
    return pixSrcX, pixSrcY




def inverseMapping(srcData, mappingImage, lens):
    I = 0
    J = 0
    mappingDict = {}
    critRad = 1.20347
    ellipticity = 0.177825
    arat = 1-ellipticity
    PA = 67.5419
    gam = 0.58
    coreRad = 0.05
    q = 0.5*critRad*pow((2-2*gam)/arat, gam)
    s = coreRad**2
    q=0.3
    # i, j are indices in image plane
    # I, J are indices in source plane
    imgXm, imgYm = commons.getGrid(-50, 50, 2, -50,  50,2, 0.05)
    deflXm = []
    deflYm = []

    for i in range(srcDimX):
        deflx, defly = fastell4py.fastelldefl(imgXm[i], imgYm[i], q, gam, arat, s)
        deflXm.append(deflx)
        deflYm.append(defly)
    new_deflXm = np.asarray(deflXm)
    new_deflYm = np.asarray(deflYm)
    srcXm = imgXm - new_deflXm
    srcYm = imgYm - new_deflYm

    plotGrid(srcXm, srcYm)



    for i in range(mappingImage.shape[1]):
        for j in range(mappingImage.shape[0]):
            I, J = getDeflectionAngle(i, j)

            mappingDict[(i,j)]= (I, J)
            if I<srcDimX and J<srcDimY:
                mappingImage[i][j]= srcData[I][J]

    return mappingDict, mappingImage


def plotGrid(xm, ym):

    #xm, ym = commons.getGrid(xSize, ySize)
    plt.plot(xm, ym, '-b')
    plt.plot(ym, xm, '-b')
    plt.show()

def main():
    dirName= "test_images/"
    imageData = commons.readFitsImage(dirName + "model_img.fits")
    srcData = commons.readFitsImage(dirName + "model_src_parameter_config.fits")
    varData = commons.readFitsImage(dirName+ "var.fits.gz")
    psfData = commons.readFitsImage(dirName+ "psf.fits.gz")
    #getModelImage()
    mappingDict, mappingImage = inverseMapping(srcData,  mappingImage=imageData, lens=0)
    commons.writeFitsImage(mappingImage, "model_test.fits")



#<<<<<<< HEAD
#Model:  SPEMD
#Critrad: 1.20347
#Ellipticity: 0.177825
#Orient_Angle: 67.5419
#Gamma: 0.58
#Core_rad: 0.05
#CenterX: 0
#CenterY: 0
#=======
#    return 0

#>>>>>>> 7d5466c18fb5bb391cfd5fc88560ddf14f779e62

if __name__=='__main__':
    main()
