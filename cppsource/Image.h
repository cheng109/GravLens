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
	long filterNPixels;
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


public:
	Image();
	Image(string imgFileName);
	void getConstants(long *filterPixelNum, long* naxis1, long* naxis2, double *res);
	void printImageInfo(int x1, int y1, int x2, int y2);
	void updateFilterImage(string regionFileName) ;
	void writeFilterImage(string imgFileName);
	void writeToFile(string imgFileName);

	virtual ~Image();
};


class Grid {
public:
	vector<double> xposList;
	vector<double> yposList;
	vector<double> briList;



};

#endif /* IMAGE_H_ */



