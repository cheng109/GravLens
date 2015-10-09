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
using namespace std;

class Image {
	string fileName;
	int naxis;
	long naxis1;
	long naxis2;
	long npixels;
	int bitpix;
	double res;
	bool filtered=0;
	//double imgXCenter;
	//double imgYCenter;

	size_t filterPixelNum;

	vector<double> data;
	vector<double> filterData;
	vector<size_t> filterX;
	vector<size_t> filterY;


public:
	Image();
	Image(string imgFileName);
	void getConstants(long *filterPixelNum, long* naxis1, long* naxis2, double *res);
	void printImageInfo(int ncol, int nrow);
	void updateFilterImage(string regionFileName) ;
	void writeFilterImage(string imgFileName);

	virtual ~Image();
};

#endif /* IMAGE_H_ */



