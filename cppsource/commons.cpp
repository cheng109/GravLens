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
		//cout << "len: " << len << endl;
		imgSize[0]=naxis1; imgSize[1] =naxis2;
		potSize[0]=naxis1; potSize[1] =naxis2;


		imgXCenter = naxis1/2.0;
		imgYCenter = naxis2/2.0;
		potXCenter = naxis1/2.0;
		potYCenter = naxis2/2.0 ;
		length = len; //dataImage->length;
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
		cout << "*********** Constants *********" << endl;
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

inline double dist(Point A, Point B) {
	return sqrt((B.x-A.x)*(B.x-A.x)+(B.y-A.y)*(B.y-A.y));
}

inline double area(Point A, Point B, Point C) {
	double side_a = dist(B, C);
	double side_b = dist(A, C);
	double side_c = dist(A, B);

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

vector<double> getTriWeight(Point A, Point B, Point C, Point P) {
	double areaA = area(P, B, C);
	double areaB = area(P, A, C);
	double areaC = area(P, A, B);
	double S = areaA + areaB + areaC;

	vector<double> w;
	w.push_back(areaA/S);
	w.push_back(areaB/S);
	w.push_back(areaC/S);
	return w;
}

double getPenalty(sp_mat* M, vec* r, vec* d, sp_mat* invC) {

	//  chi2 =  (M*r-d)^T(M*r-d)
	//cout << *M << endl;
	vec res = (*M)*(*r)-(*d);
	vec chi2 =  res.t()*(*invC)*res;

	return chi2(0,0);
}


void getLinearInterpolate(Point A, Point B,  Point C,  Point *P,  char direction) {

	double a=0, b=0; 

	if(direction=='x') {
		if(abs(A.x-B.x)<10e-8)  {
			P->x = 0.5*(A.x+B.x);
			P->y = C.y;
		}
		else {
			a = (A.y-B.y)/(A.x-B.x);
			b = A.y-a*A.x;
			P->x = (C.y-b)/a;
			P->y = C.y;
		}
	}
	if(direction=='y') {

		a = (A.y-B.y)/(A.x-B.x);
		b = A.y-a*A.x;
		P->x = C.x;
		P->y = a*C.x+b;
	}

}



vector<double> getPentWeigth(Point A, Point B, Point C, Point D, Point E) {
	vector<double> pentWeight;

	Point P(0, 0, 0);
	Point Q(0, 0, 0);
	Point M(0, 0, 0);
	Point N(0, 0, 0);
	getLinearInterpolate(A, B, C, &Q, 'x');
	getLinearInterpolate(D, E, C, &P, 'x');

	getLinearInterpolate(B, E, C, &M, 'y');
	getLinearInterpolate(A, D, C, &N, 'y');

	double dCQ_AB = dist(C, Q)*dist(A, B);
	cout << "Dist(D, C) " << Q.x << "\t" << Q.y  << endl;
	double dCP_DE = dist(C, P)*dist(D, E);
	//double dAB = ;

	pentWeight.push_back(dist(Q, B)/dCQ_AB);    // wAx
	pentWeight.push_back(dist(Q, A)/dCQ_AB); 		// wBx
	pentWeight.push_back(-(1/dist(C, P)+1/dist(C, Q)));  //wCx
	pentWeight.push_back(dist(P, E)/dCP_DE);   //wDx
	pentWeight.push_back(dist(P, D)/dCP_DE); 			//wEx

	double dCN_AD = dist(C, N)*dist(A, D);
	double dCM_BE = dist(C, M)*dist(B, E);

	pentWeight.push_back(dist(N, D)/dCN_AD);   //wAy
	pentWeight.push_back(dist(M, E)/dCM_BE);   //wBy
	pentWeight.push_back(-(1/dist(C, N)+1/dist(C,M)));  //wCy
	pentWeight.push_back(dist(A, N)/dCN_AD);    //wDy
	pentWeight.push_back(dist(B, M)/dCM_BE);    // wEy

	return pentWeight;

}



normVec getNormVector(Point A, Point B, Point C) {
	normVec norm; 
	norm.n0 = (C.y-B.y)*(B.z-A.z)-(C.z-B.z)*(B.y-A.y); 
	norm.n1 = (C.z-B.z)*(B.x-A.x)-(C.x-B.x)*(B.z-A.z); 
	norm.n2 = (C.x-B.x)*(B.y-A.y)-(C.y-B.y)*(B.x-A.x); 
	double s = sqrt(norm.n0*norm.n0 + norm.n1*norm.n1 + norm.n2*norm.n2);
	norm.n0 = norm.n0/s;
	norm.n1 = norm.n1/s;
	norm.n2 = norm.n2/s;

	return norm;
}

normVec meanNormVector(vector<normVec>  normList) {
	normVec meanNorm(0, 0, 0);
	for(int i=0; i<normList.size(); ++i) {
		meanNorm.n0 += normList[i].n0;
		meanNorm.n1 += normList[i].n1;
		meanNorm.n2 += normList[i].n2;
	}

	double s = meanNorm.n0*meanNorm.n0
			+ meanNorm.n1*meanNorm.n1
			+ meanNorm.n2*meanNorm.n2;

	if(s==0) {
		meanNorm.n0 = 0;
		meanNorm.n1 = 0;
		meanNorm.n2 = 0;
	}
	else {
		meanNorm.n0 = meanNorm.n0/s;
		meanNorm.n1 = meanNorm.n1/s;
		meanNorm.n2 = meanNorm.n2/s;
	}
	return meanNorm;

}





