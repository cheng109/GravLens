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
#include <iomanip>
#include <armadillo>
#include <cmath>
using namespace arma;
using namespace std;
#define buffsize 10000


Image::Image() {
	// TODO Auto-generated constructor stub

}

Image::Image(vector<double> xpos, vector<double> ypos, vector<double> *briList, long naxis1, long naxis2, int bitpix):
			naxis(2), naxis1(naxis1), naxis2(naxis2), bitpix(bitpix), data(naxis1*naxis2, 0) {

	npixels = naxis1*naxis2;
	long iList=0, x, y;
	for(int i=0; i< briList->size(); ++i) {
		x = nearbyint(xpos[i]);
		y = nearbyint(ypos[i]);

		if(x>0 && x< naxis1 && y>0 && y<naxis2) {
			iList = naxis1*y+x;
		}
		//cout <<i << "\t" << xpos[i] << "\t" << ypos[i] << "\t" << iList << "\t" << endl;

		data[iList] += briList->at(i);
			
	}

}


Image::Image(string imgFileName) {

	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	int status,  nfound, anynull;
	long naxes[2], fpixel, nbuffer;


	float nullval, buffer[buffsize];
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

	while (npixels > 0)
	{
		nbuffer = npixels;
		if (npixels > buffsize)
			nbuffer = buffsize;     /* read as many pixels as will fit in buffer */

		if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval,
				buffer, &anynull, &status) )
			printerror( status );
		for(long i=0; i<nbuffer; ++i) {
			data.push_back(buffer[i]);
		}
		npixels -= nbuffer;    /* increment remaining number of pixels */
		fpixel  += nbuffer;    /* next pixel to be read in image */
	}

	npixels  = naxes[0] * naxes[1];
	if ( fits_close_file(fptr, &status) )
		printerror( status );
}



void Image::getConstants(long *len, long* naxis1, long* naxis2, double *res, int *bit){
	cout << "length " << length << endl;
	*len=this->length;
	*naxis1 = this->naxis1;
	*naxis2 = this->naxis2;
	*res = this->res;
	*bit = this->bitpix;
}



void Image::printImageInfo(int x1, int y1, int x2, int y2) {
	cout << "************************" << endl;
	cout << "File:   " << fileName << endl;
	cout << "NAXIS:  " << naxis << endl;
	cout << "NAXIS1: " << naxis1 << endl;
	cout << "NAXIS2: " << naxis2 << endl;
	cout << "Total:  " << npixels << endl;
	cout << "TYPE:   " << bitpix << endl;
	cout << "************************" << endl;

	//cout << counter << endl;
	x1 = (x1>0)?x1:0;
	y1 = (y1>0)?y1:0;
	x2 = (x2<naxis1)?x2:naxis1;
	y2 = (y2<naxis2)?y2:naxis2;
	// in order of DS9 show.
	for (int y=y2; y>y1; --y) {
		for(int x=x1; x<x2; ++x) {
			cout <<setprecision(4) << this->data[y*naxis1+x] << "\t";
		}
		cout << endl;

	}
	cout << endl;
	cout << data.size() << endl;
}


void Image::writeFilterImage(string imgFileName) {
	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	int status, ii;
	long  fpixel;
	long naxes[2] = { naxis1, naxis2 };   /* image is 300 pixels wide by 200 rows */
	remove(imgFileName.c_str());               /* Delete old file if it already exists */
	status = 0;         /* initialize status before calling fitsio routines */
	if (fits_create_file(&fptr, imgFileName.c_str(), &status)) /* create new FITS file */
		printerror( status );           /* call printerror if error occurs */
	if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
		printerror( status );
	double **array = (double **)malloc(naxis2* sizeof(double**));
	array[0] = (double *)malloc( naxis1 * naxis2* sizeof( double ) );
	for( ii=1; ii<naxis2; ii++ )
		array[ii] = array[ii-1] + naxis1;
	for (int i = 0; i < naxis2; i++) {
		for (int j = 0; j < naxis1; j++)
			array[i][j] = 0.0;
	}
	for( ii=1; ii<naxis2; ii++ )
		array[ii] = array[ii-1] + naxis1;
	cout << dataList.size() << endl;
	for (int i=0; i<dataList.size();  ++i) {
		array[yList[i]][xList[i]] = dataList[i];
	}
	fpixel = 1;                               /* first pixel to write      */
	if (fits_write_img(fptr, TDOUBLE, fpixel, npixels, array[0], &status))
		printerror( status );
	free( array[0] );  /* free previously allocated memory */
	free( array);
	if ( fits_close_file(fptr, &status) )                /* close the file */
		printerror( status );


}


void Image::writeToFile(string imgFileName) {
	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	int status, ii;
	long  fpixel;

	long naxes[2] = { naxis1, naxis2 };   /* image is 300 pixels wide by 200 rows */
	remove(imgFileName.c_str());               /* Delete old file if it already exists */
	status = 0;         /* initialize status before calling fitsio routines */
	if (fits_create_file(&fptr, imgFileName.c_str(), &status)) /* create new FITS file */
		printerror( status );           /* call printerror if error occurs */
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

void Image::updateFilterImage(string regionFileName) {
	vector<double> xpos, ypos;
	int lenRegionList =parseReagionFile(regionFileName, &xpos, &ypos);

	for(int x=0; x<naxis1; ++x) {
		for (int y=0; y<naxis2; ++y) {
		  if(pnpoly(lenRegionList, &xpos, &ypos,  double(x+0.5) ,  double(y+0.5))) {

				dataList.push_back(data[naxis1*y+x]);
				xList.push_back(x+1);
				yList.push_back(y+1);
			}
		}
	}

	length = dataList.size() ;
	for(int i=0; i<dataList.size(); ++i) {
		iList.push_back(i);
	}
	//cout << "length" << "\t" << length<<endl;
}


vec Image::getMatrixD() {
	vec d(length);
	for(int i=0; i<length; ++i) {
		d[i]= dataList[i];
	}
	return d;

}
void Image::updateGridPointType() {
	for(int i=0; i<length; ++i) {
		if((xList[i]+yList[i])%2==0)
			type.push_back(0);
		else
			type.push_back(1);
	}

}

sp_mat Image::getVarMatrix()  {

	//this->updateFilterImage(regionFileName);
	long n= length;
	//n = 100;
	sp_mat invC = speye<sp_mat>(n, n);
	for(int i=0; i<n; ++i) {
		invC(i,i)=varList[i];
	}
	return invC;
}

void Image::updateVarList(double threshold, double backVar) {
	for(int i=0; i<length; ++i) {
		if(dataList[i]<threshold) {
			varList.push_back(backVar);
		}
		else
			varList.push_back(dataList[i]);
	}
}

void Image::updateVarList(string varFileName, string regionFileName) {
	Image varImage(varFileName);
	varImage.updateFilterImage(regionFileName);
	for(int i=0; i<length; ++i) {
		varList.push_back(dataList[i]);
	}
}

/*

void Image::erasePixel(int index) {

	dataList.erase(dataList.begin()+index);
	varList.erase(varList.begin()+index);
	xList.erase(xList.begin()+index);
	yList.erase(yList.begin()+index);

	type.erase(type.begin()+index);
	iList.erase(iList.begin()+index);
	length = length -1;

}
*/


Image::~Image() {
	// TODO Auto-generated destructor stub
}
