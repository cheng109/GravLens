__author__ = 'juncheng'
import commons
from fastell4py import fastell4py
import math
import models

def getDeflection(imgX, imgY, model, critRad, const):
    imgX = (imgX-const.imgXCenter)*const.imgRes
    imgY = (imgY-const.imgYCenter)*const.imgRes

    if model=='SPEMD':
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


    if model=='SIE':
        # rotation by angles
        core=0.001
        x1 = imgX
        y1 = imgY
        #critRad = 1.25
        LM_SIE = models.SIE(critRad=critRad,axisRatio=1, pa=0, centerX=0, centerY=0)
        q = LM_SIE.axisRatio
        if q==1:
            q=0.999
        root1mq = math.sqrt(1.0-q**2)
        phi = math.sqrt(q**2*(core**2+x1**2)+y1**2)
        fac = LM_SIE.critRad*math.sqrt(q)/root1mq
        deflx = fac*math.atan(root1mq*x1/(phi+core))
        defly = fac*commons.lm_arctanh(root1mq*y1/(phi+core*LM_SIE.axisRatio**2))
        # Rotate back
        pdeflx = deflx
        pdefly = defly


    srcX = (imgX - pdeflx)/const.srcRes+const.srcXCenter
    srcY = (imgY - pdefly)/const.srcRes+const.srcYCenter
    return deflx, defly, srcX, srcY






if __name__=='__main__':
    main()
