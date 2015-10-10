/*
 * Image.h
 *
 *  Created on: Oct 8, 2015
 *      Author: cheng109
 */

#ifndef IMAGE_H_
#define IMAGE_H_


#include<string>
#include<vector>
#include <armadillo>

using namespace std;
using namespace arma;

class Const;

class Image {
friend class Model;
	string fileName;
	int naxis;
	long naxis1;
	long naxis2;
	long npixels;

	int bitpix;
	double res;
	//bool filtered=0;
	//double imgXCenter;
	//double imgYCenter;

	size_t filterPixelNum;

	vector<double> data;
	vector<double> filterData;
	vector<int> filterX;
	vector<int> filterY;
	vector<int> index;


public:
	Image();
	Image(vector<double> xpos, vector<double> ypos,vector<double> *briList, long naxis1, long naxis2, int bitpix);

	Image(string imgFileName);
	void getConstants(long *filterPixelNum, long* naxis1, long* naxis2, double *res, int* bit);
	void printImageInfo(int x1, int y1, int x2, int y2);
	void updateFilterImage(string regionFileName) ;
	void writeFilterImage(string imgFileName);
	void writeToFile(string imgFileName);

	vec getMatrixD();

	virtual ~Image();
};


#endif /* IMAGE_H_ */



