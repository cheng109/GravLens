/*
 * commons.cpp
 *
 *  Created on: Sep 30, 2015
 *      Author: cheng109
 */
#include "commons.h"
#include "fitsio.h"
#include <string>
#include "Image.h"
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <armadillo>
#include <map>


#define buffsize 1000
using namespace arma;
using namespace std;



Conf::Conf(Image* dataImage, map<string, string> confMap) {
		//imgSize[0] = dataImage->
		long naxis1, naxis2, len ;
		int bit;
		double res;
		dataImage->getConstants(&len, &naxis1, &naxis2, &res, &bit);
		res = 0.3;

		imgSize[0]=naxis1; imgSize[1] =naxis2;
		potSize[0]=naxis1; potSize[1] =naxis2;


		imgXCenter = naxis1/2.0;
		imgYCenter = naxis2/2.0;
		potXCenter = naxis1/2.0;
		potYCenter = naxis2/2.0 ;
		length = len;
		bitpix = bit;

		srcRes = stod(confMap["srcRes"]);
		imgRes = stod(confMap["imgRes"]);
		potRes = stod(confMap["potRes"]);

		srcSize[0] =stod(confMap["srcX"]);
		srcSize[1] =stod(confMap["srcY"]);

		srcXCenter = srcSize[0]/2.0;
		srcYCenter = srcSize[1]/2.0;


}



void Conf::printConfList(){
		cout << "*********** ConfANTS *********" << endl;
		cout << "srcSize:    " << srcSize[0] << ",\t"<< srcSize[1] << endl;
		cout << "imgSize:    " << imgSize[0] << ",\t"<<imgSize[1]  << endl;
		cout << "potSize:    " << potSize[0] << ",\t"<<potSize[1]  << endl;
		cout << "srcRes:     " << srcRes << endl;
		cout << "imgRes:     " << imgRes << endl;
		cout << "potRes:     " << potRes << endl;
		cout << "srcCenter:  " << srcXCenter << ",\t"<<srcYCenter << endl;
		cout << "imgCenter:  " << imgXCenter << ",\t"<<imgYCenter << endl;
		cout << "potCenter:  " << potXCenter << ",\t"<<potYCenter << endl;
		cout << "length:     " << length << endl;
		cout << "*******************************" << endl;
}

void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/
    if (status)
    {
       fits_report_error(stderr, status); /* print error report */
       exit( status );    /* terminate the program, returning error status */
    }
    return;
}

bool pnpoly(size_t nvert, vector<double> *vertx, vector<double> *verty, double testx, double testy)
{
  int i, j;
  bool c = FALSE;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty->at(i)>testy) != (verty->at(j)>testy)) &&
     (testx < (vertx->at(j)-vertx->at(i)) * (testy-verty->at(i)) / (verty->at(j)-verty->at(i)) + vertx->at(i)) )
       c = !c;
  }
  return c;
}


int parseReagionFile(string regionFileName, vector<double> *xpos, vector<double> *ypos) {

	ifstream regionFile(regionFileName.c_str());
	string line, token;
	size_t pos=0;
	//vector<double> xpos;
	//vector<double> ypos;

	while (getline(regionFile, line)) {
		if (line[0]!='#' && line.substr(0, 6)!="global" && line.substr(0, 5)!="image") {
			size_t pos1 = line.find_first_of("(", pos);
			size_t pos2 = line.find_first_of(")", pos);

			istringstream ss(line.substr(pos1+1, pos2-pos1));
			int flag=-1;
			while(getline(ss, token, ',' )){
				if(flag<0) xpos->push_back(stod(token));
				if(flag>0) ypos->push_back(stod(token));
				//cout << token << endl;
				flag = (-1)*flag;
			};
		}
		if(xpos->size()!=ypos->size())
			cout << "Error when reading region File!" << endl;
	}
	return xpos->size();

}


map<string, string> parseConfigure(string confFileName) {

	map<string, string> confMap;
	ifstream confFile(confFileName.c_str());
	string line;
	vector<string> items;
	while (getline(confFile, line)) {
		//cout << line << endl;
		items.clear();
		if(line[0]!='#') {
			istringstream iss(line);
			//cout << iss << endl;
			copy(istream_iterator<string>(iss),
			      istream_iterator<string>(),
			      back_inserter(items));
			if (items.size()>=2)	  {
				confMap[items[0]] = items[1];
				cout << items[0] << "\t" << items[1] << endl;
			}
		}
	}

	return confMap;

}


inline double dist(double Ax, double Ay, double Bx, double By) {
	return sqrt((Bx-Ax)*(Bx-Ax)+(By-Ay)*(By-Ay));
}

inline double area(double Ax, double Ay, double Bx, double By, double Cx, double Cy) {
	double side_a = dist(Bx, By, Cx, Cy);
	double side_b = dist(Ax, Ay, Cx, Cy);
	double side_c = dist(Ax, Ay, Bx, By);

	double s = 0.5*(side_a + side_b + side_c);
	return sqrt(s*(s-side_a)*(s-side_b)*(s-side_c));
}


double lm_arctanh(double x) {
	if (x < -1 || x > 1.) {
		fprintf(stderr,"lm_arctanh: invalid x: %g. Must be 0 <= x <= 1\n",x);
		return 0;
	}
	return log(sqrt((1.+x)/(1.-x)));
}



vector<double> getTriWeight(double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Px, double Py) {
	double areaA = area(Px, Py, Bx, By, Cx, Cy);
	double areaB = area(Px, Py, Ax, Ay, Cx, Cy);
	double areaC = area(Px, Py, Ax, Ay, Bx, By);
	double S = areaA + areaB + areaC;

	vector<double> w;
	w.push_back(areaA/S);
	w.push_back(areaB/S);
	w.push_back(areaC/S);

	return w;


}


double getPenalty(sp_mat* M, vec* r, vec* d, sp_mat* C) {

	//  chi2 =  (M*r-d)^T(M*r-d)
	//cout << *M << endl;
	vec res = (*M)*(*r)-(*d);
	vec chi2 =  res.t()*res;

	return chi2(0,0);

}







