__author__ = 'juncheng'
import commons
from fastell4py import fastell4py


def getDeflection(imgX, imgY, model,const):

    if model=='SPEMD':
        critRad = 1.20347*10
        ellipticity = 0.177825
        arat = 1-ellipticity
        PA = 67.5419
        gam = 0.58
        coreRad = 0.05
        q = 0.3#0.5*critRad*pow((2-2*gam)/arat, gam)
        s = coreRad**2

        imgX = (imgX-const.imgXCenter)*const.imgRes
        imgY = (imgY-const.imgYCenter)*const.imgRes

       # angleImageX = (pixImageX-61)*imgRes
       # angleImageY = (pixImageY-61)*imgRes

#        pixSrcX=int( (angleImageX - deflx)/srcRes + srcDimX/2)
 #       pixSrcY= int ((angleImageY - defly)/srcRes + srcDimY/2)

        deflx, defly= fastell4py.fastelldefl(imgX, imgY,  q, gam, arat, s)
        srcX = (imgX - deflx)/const.srcRes+const.srcXCenter
        srcY = (imgY - defly)/const.srcRes+const.srcYCenter


    return deflx, defly, srcX, srcY






if __name__=='__main__':
    main()
