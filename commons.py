
from astropy.io import fits


def readFitsImage(imageName):
    hdulist = fits.open(imageName)
    data =  hdulist[0].data
    hdulist.close()
    return data

def writeFitsImage(data, outputName):
    hdu = fits.PrimaryHDU(data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(outputName, clobber=True)


if __name__=='__main__':
    main()
