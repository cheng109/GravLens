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
#include <map>
#include <armadillo>

using namespace std;
using namespace arma;

class Conf{
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
	long potN;
	int bitpix;

	Conf(Image* image, map<string, string> confMap);
	void printConfList();
};


inline double dist(double Ax, double Ay, double Bx, double By) ;
inline double area(double Ax, double Ay, double Bx, double By, double Cx, double Cy);
vector<double> getTriWeight(double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Px, double Py);

void printerror( int status);

bool pnpoly(size_t nvert, vector<double> *vertx, vector<double> *verty, double testx, double testy);
void updateConf(string confFileName);
int parseReagionFile(string regionFileName, vector<double> *xpos, vector<double> *ypos);
map<string, string> parseConfigure(string confFileName) ;
double getPenalty(sp_mat* M, vec* r, vec* d, sp_mat* C);
double lm_arctanh(double x);

#endif
