/*
 * Model.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: juncheng
 */

#include "Model.h"
#include <armadillo>
#include <iostream>
#include <string>
#include <map>
#include <cmath>
using namespace std;
using namespace arma;

Model::Model() {
	// TODO Auto-generated constructor stub

}


Model::Model(string name, double modelCenterX, double modelCenterY, double critR, double e,  double PA, double mass):
		name(name), modelCenterX(modelCenterX), modelCenterY(modelCenterY), critR(critR), e(e), PA(PA), mass(mass) {



}



void Model::getDeflectionAngle(Const* conList, int imgX, int imgY, double *srcX, double *srcY) {
	double fX = (imgX -conList->imgXCenter-modelCenterX)*conList->imgRes;   // Model center frame
	double fY = (imgY -conList->imgYCenter-modelCenterY)*conList->imgRes;

	double pfX = (imgX-conList->imgXCenter)*conList->imgRes; 			// image center frame;
	double pfY = (imgY-conList->imgYCenter)*conList->imgRes;

	double pDeltaX = 0;
	double pDeltaY = 0;

	double fDenom;

	if(name.compare("PTMASS")==0) {
		fDenom = fX*fX+fY*fY;
		double fMult = critR*critR/fDenom;

		pDeltaX = fX*fMult;
		pDeltaY = fY*fMult;

	}

	if(name.compare("SIE")==0) {
		cout << endl;
	}


	*srcX = (pfX - pDeltaX)/conList->srcRes+conList->srcXCenter;
	*srcY = (pfY - pDeltaY)/conList->srcRes+conList->srcYCenter;
	if (*srcX >70)
	cout << imgX << "\t" << imgY << "\t" << fX << "\t" << fY << "\t"<<sqrt(fDenom) << "\t"<< pDeltaX <<"\t" << pDeltaY <<"\t" << *srcX << "\t" << *srcY <<endl;

}
Model::~Model() {
	// TODO Auto-generated destructor stub
}


map<pair<int, int>,int> Model::createPosMapping(Image* image, vector<double>* srcXList,vector<double>* srcYList, Const* conList) {

	map<pair<int, int>,int> posMap;

	for(int i=0; i<conList->length; ++i) {
		int imgX = image->filterX[i];
		int imgY = image->filterY[i];
		getDeflectionAngle(conList,imgX, imgY,  &srcXList->at(i),  &srcYList->at(i));

		//if (srcXList->at(i)>60)
		//cout << srcXList->at(i) << "\t" << srcYList->at(i) << endl;
		posMap[make_pair(imgX, imgY)] = i;
	}
	return posMap;
}

mat buildLensMatrix(Image* image, Model* model, Const* constList) {

	mat L(constList->length,constList->length, fill::eye);

	return L;

}
