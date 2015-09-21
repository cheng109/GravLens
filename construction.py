__author__ = 'cheng109'

import numpy as np
import commons
## construction for source image vector



def getLensOperator(mappingDict):
    dim = len(mappingDict)
    L = np.zeros((dim, dim))
    imgPointList = mappingDict.keys()
    srcPointList = mappingDict.values()
    srcPosition = []
    d =  np.zeros(dim)
    for i in range(dim):
        imgX, imgY = imgPointList[i]
        srcX, srcY, imgBrightNess, type = srcPointList[i]
        srcPosition.append((srcX, srcY))
        d[i] = imgBrightNess

        if type=='o' and (imgX, imgY-1) in imgPointList and (imgX, imgY+1) in imgPointList and (imgX+1, imgY) in imgPointList:
            w1Index = imgPointList.index((imgX, imgY-1))
            w2Index = imgPointList.index((imgX, imgY+1))
            w3Index = imgPointList.index((imgX+1, imgY))
            A = (mappingDict[(imgX,imgY-1)][0], mappingDict[(imgX,imgY-1)][1])
            B = (mappingDict[(imgX,imgY+1)][0], mappingDict[(imgX,imgY+1)][1])
            C = (mappingDict[(imgX+1,imgY)][0], mappingDict[(imgX+1,imgY)][1])

            L[i][w1Index], L[i][w2Index],L[i][w3Index] = commons.getTriWeight(A,B,C, (srcX, srcY))
            #print L[i][w1Index], L[i][w2Index],L[i][w3Index]
        else:
            L[i][i] =1

    L = np.identity(dim)
    s = np.linalg.solve(L,d)
    # Now we have the 'SrcPosition' and 'SrcBrightness'

    return srcPosition, s



def main():



    return


if __name__ =="__main__":
    main()