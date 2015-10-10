/*
 * commons.h
 *
 *  Created on: Sep 30, 2015
 *      Author: cheng109
 */

#ifndef COMMONS_H_
#define COMMONS_H_


#include <string>
#include <vector>
#include "Image.h"
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
	double potXCenter;
	double potYCenter;
	double length;
	int bitpix;

	Const(Image* image);
	void printConstList();
};


void printerror( int status);

bool pnpoly(size_t nvert, vector<double> *vertx, vector<double> *verty, double testx, double testy);

int parseReagionFile(string regionFileName, vector<double> *xpos, vector<double> *ypos);


#endif
