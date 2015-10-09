/*
 * Image.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: cheng109
 */

#include "Image.h"
#include "fitsio.h"
#include <string>
#include <vector>
#include <iostream>
#include "commons.h"
using namespace std;
#define buffsize 10000


Image::Image() {
	// TODO Auto-generated constructor stub

}

Image::Image(string imgFileName) {

	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	int status,  nfound, anynull;
	long naxes[2], fpixel, nbuffer,  ii;


	float datamin, datamax, nullval, buffer[buffsize];
	//char filename[]  = "atestfil.fit";     /* name of existing FITS file   */

	status = 0;


	if ( fits_open_file(&fptr, imgFileName.c_str(), READONLY, &status) )
		printerror( status );
	if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
		printerror( status );
	if(fits_get_img_type(fptr, &bitpix, &status))
		printerror( status );
	if (fits_get_img_dim(fptr, &naxis, &status))
		printerror( status );
	naxis1 = naxes[0] ;
	naxis2 = naxes[1] ;
	npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
	fpixel   = 1;
	nullval  = 0;                /* don't check for null values in the image */
	datamin  = 1.0E30;
	datamax  = -1.0E30;

	while (npixels > 0)
	{
		nbuffer = npixels;
		if (npixels > buffsize)
			nbuffer = buffsize;     /* read as many pixels as will fit in buffer */

		if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval,
				buffer, &anynull, &status) )
			printerror( status );
		for(long i=0; i<nbuffer; ++i)
			data.push_back(buffer[i]);

		cout << endl;
		npixels -= nbuffer;    /* increment remaining number of pixels */
		fpixel  += nbuffer;    /* next pixel to be read in image */
	}

	npixels  = naxes[0] * naxes[1];
	if ( fits_close_file(fptr, &status) )
		printerror( status );
}


void Image::printImageInfo(int ncol, int nrow) {
	cout << "************************" << endl;
	cout << "File:   " << fileName << endl;
	cout << "NAXIS:  " << naxis << endl;
	cout << "NAXIS1: " << naxis1 << endl;
	cout << "NAXIS2: " << naxis2 << endl;
	cout << "Total:  " << npixels << endl;
	cout << "TYPE:   " << bitpix << endl;
	cout << "************************" << endl;

	//cout << counter << endl;
	for(int i=0; i<npixels; ++i){
		cout << data.at(i) << "\t";
	}
	cout << endl;
}


void Image::writeToFile(string imgFileName) {
	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	int status, ii, jj;
	long  fpixel;

	/* initialize FITS image parameters */

	//char filename[] = "";             /* name for new FITS file */
	long naxes[2] = { naxis1, naxis2 };   /* image is 300 pixels wide by 200 rows */


	remove(imgFileName.c_str());               /* Delete old file if it already exists */

	status = 0;         /* initialize status before calling fitsio routines */

	if (fits_create_file(&fptr, imgFileName.c_str(), &status)) /* create new FITS file */
		printerror( status );           /* call printerror if error occurs */

	/* write the required keywords for the primary array image.     */
	/* Since bitpix = USHORT_IMG, this will cause cfitsio to create */
	/* a FITS image with BITPIX = 16 (signed short integers) with   */
	/* BSCALE = 1.0 and BZERO = 32768.  This is the convention that */
	/* FITS uses to store unsigned integers.  Note that the BSCALE  */
	/* and BZERO keywords will be automatically written by cfitsio  */
	/* in this case.                                                */

	if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
		printerror( status );

	double **array = (double **)malloc(naxis2* sizeof(double**));
	array[0] = (double *)malloc( naxis1 * naxis2* sizeof( double ) );
	for( ii=1; ii<naxis2; ii++ )
	      array[ii] = array[ii-1] + naxis1;
	for (int i = 0; i < naxis2; i++) {
		for (int j = 0; j < naxis1; j++){
			array[i][j] = data[i*naxis1+j];
		}
	}
	fpixel = 1;                               /* first pixel to write      */
	/* write the array of unsigned integers to the FITS file */
	if (fits_write_img(fptr, TDOUBLE, fpixel, npixels, array[0], &status))
		printerror( status );

	free( array[0] );  /* free previously allocated memory */
	free( array);
	if ( fits_close_file(fptr, &status) )                /* close the file */
		printerror( status );
}


Image::~Image() {
	// TODO Auto-generated destructor stub
}

