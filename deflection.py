__author__ = 'juncheng'
import commons
import numpy as np
from fastell4py import fastell4py
import math
import models

def getDeflection(imgX, imgY, MODEL, const):
     # coordination center at (lensMid, lensMid)
    fX = (imgX-const.imgXCenter-MODEL.centerX)*const.imgRes
    fY = (imgY-const.imgYCenter-MODEL.centerY)*const.imgRes

    # coordination center at (imgMid, imgMid)
    pfX = (imgX-const.imgXCenter)*const.imgRes
    pfY = (imgY-const.imgYCenter)*const.imgRes

    if MODEL.name=='PTMASS':
        fDenom = fX**2+fY**2
        fMult = MODEL.critRad**2/fDenom

        pDeltaX = fX*fMult
        pDeltaY = fY*fMult


    if MODEL.name=='SIE':
        # rotation by angles

        fCore=0.0
        MODEL.pa = np.radians(MODEL.pa)
        fq = MODEL.axisRatio
        fCosTheta = np.cos(MODEL.pa)
        fSinTheta = np.sin(MODEL.pa)
        if fq==1.0:
            fq=0.999
        x1 = fX*fCosTheta + fY*fSinTheta
        y1 = -fX*fSinTheta + fY*fCosTheta

        root1mq = math.sqrt(1.0-fq*fq)

        phi = np.sqrt(fq*fq*(fCore*fCore+x1*x1)+y1*y1)   #original

        fac = MODEL.critRad*np.sqrt(fq)/root1mq

        deltax1 = fac*math.atan(root1mq*x1/(phi+fCore))
        deltay1 = fac*commons.lm_arctanh(root1mq*y1/(phi+fCore*fq*fq))
        # Rotate back
        pDeltaX = deltax1 * fCosTheta - deltay1 * fSinTheta
        pDeltaY = deltay1 * fCosTheta + deltax1 * fSinTheta


    if MODEL.name=='SPEMD':
        critRad = 1.20347*10
        ellipticity = 0.177825
        arat = 1-ellipticity
        PA = 67.5419
        gam = 0.58
        coreRad = 0.05
        q = 50.03#0.5*critRad*pow((2-2*gam)/arat, gam)
        s = coreRad**2
        deflx, defly= fastell4py.fastelldefl(imgX, imgY,  q, gam, arat, s)
        pdeflx = deflx
        pdefly = defly


    srcX = (pfX - pDeltaX)/const.srcRes+const.srcXCenter
    srcY = (pfY - pDeltaY)/const.srcRes+const.srcYCenter


    #print "imgX= ", imgX, "imgY= ", imgY, "pfX= ", pfX, "pfY= ", pfY, "fX= ", fX, "fY= ", fY, "R= ", np.sqrt(fX**2+fY**2), "pDeltaX = ", pDeltaX, "pDeltaY =", pDeltaY, "srcX = " , srcX, "srcY = ", srcY
    return  pDeltaX, pDeltaY, srcX, srcY






if __name__=='__main__':
    main()
