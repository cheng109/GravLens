/*
 * commons.h
 *
 *  Created on: Sep 30, 2015
 *      Author: cheng109
 */

#include <string>
#include <vector>
#include "Image.h""
using namespace std;

class Const{
public:
	size_t srcSize[2];
	size_t imgSize[2];
	size_t potSize[2];
	double srcRes;
	double imgRes;
	double potRes;
	double srcXCenter;
	double srcYCenter;
	double imgXCenter;
	double imgYCenter;
	double length;

	Const(Image* dataImage) {
		//imgSize[0] = dataImage->
		dataImage->getConstants(length, long* naxis1, long* naxis2, double *res);
		srcSize[0]=1; srcSize[1] =1;
		imgSize[0]=1; imgSize[1] =2;
		potSize[0]=1; potSize[1]=2;
		srcRes =  dataImage->res;
		imgRes = dataImage->res;
		potRes = dataImage->res;
		srcXCenter = 1;
		srcYCenter =2;
		imgXCenter = 1;
		imgYCenter =2 ;
		length = dataImage->npixels;


	}
};



void printerror( int status);

bool pnpoly(size_t nvert, vector<double> *vertx, vector<double> *verty, double testx, double testy);

size_t parseReagionFile(string regionFileName, vector<double> *xpos, vector<double> *ypos);
