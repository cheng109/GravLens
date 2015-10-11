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



void Model::getDeflectionAngle(Conf* conList, int imgX, int imgY, double *srcX, double *srcY) {
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
		//if (fX == 0 && fY == 0) {
		//	*pDeltaX = *pDeltaY = pLensComp->fParameter[0];
		//
		//						break;
		//					}
/*		real_t	phi,root1mq,fq,fac,fCore=0,fCosTheta,fSinTheta,x1,y1,deltax1,deltay1;



						/* pre-calculate constants
						fCosTheta = cos(pLensComp->fParameter[2]*M_PI/180);
						fSinTheta = sin(pLensComp->fParameter[2]*M_PI/180);

						fq = pLensComp->fParameter[1];
						  if (fq > 1.0) {
						                    iStatus = LM_BADPROJ;
						                    break;
						                }
						                if (fq==1.0) fq = 0.999;

										 rotate reference frame to x-axis
										x1 = fX*fCosTheta + fY*fSinTheta;
										y1 = -fX*fSinTheta + fY*fCosTheta;

										root1mq = sqrt(1.0-fq*fq);
										phi = sqrt(fq*fq*(fCore*fCore + x1*x1) + y1*y1);
										 use sqrt(fq) here not fq. This presevers scale of Einstein Radius
										fac = pLensComp->fParameter[0]*sqrt(fq)/root1mq;

										deltax1 = fac*atan(root1mq*x1/(phi + fCore));
										deltay1 = fac*lm_arctanh(root1mq*y1/(phi+ fCore*fq*fq));

										 rotate back again
										*pDeltaX = deltax1*fCosTheta - deltay1*fSinTheta;
										*pDeltaY = deltay1*fCosTheta + deltax1*fSinTheta;

						                 printf("x1,y1: %g,%g. root1mq: %g, phi: %g, fac: %g\n",x1,y1,root1mq,phi,fac);

									}
				*/
	}


	*srcX = (pfX - pDeltaX)/conList->srcRes+conList->srcXCenter;
	*srcY = (pfY - pDeltaY)/conList->srcRes+conList->srcYCenter;
	//	if (*srcX >70)
	//cout << imgX << "\t" << imgY << "\t" << fX << "\t" << fY << "\t"<<sqrt(fDenom) << "\t"<< pDeltaX <<"\t" << pDeltaY <<"\t" << *srcX << "\t" << *srcY <<endl;

}
Model::~Model() {
	// TODO Auto-generated destructor stub
}


map<pair<int, int>,int> Model::createPosMapping(Image* image, vector<double>* srcXList,vector<double>* srcYList, Conf* conList) {

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

mat buildLensMatrix(Image* image, Model* model, Conf* constList) {

	mat L(constList->length,constList->length, fill::eye);

	return L;

}
