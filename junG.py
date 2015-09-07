
import commons


def getModelImage(src, lens):

    model = 0
    return model


def main():
    dirName= "test_images/"

    imageData = commons.readFitsImage(dirName + "model_img.fits")
    srcData = commons.readFitsImage(dirName + "model_src.fits")
    varData = commons.readFitsImage(dirName+ "var.fits.gz")
    psfData = commons.readFitsImage(dirName+ "psf.fits.gz")



    return 0



if __name__=='__main__':
    main()
