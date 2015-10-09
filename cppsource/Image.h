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

	vector<double> data;

	long counter;



public:
	Image();
	Image(string imgFileName);
	void writeToFile(string imgFileName);
	void printImageInfo(int ncol, int nrow);





	virtual ~Image();
};

#endif /* IMAGE_H_ */
