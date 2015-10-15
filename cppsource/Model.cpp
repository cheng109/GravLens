/*
 * Model.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: juncheng
 */

#include "Model.h"
#include "commons.h"
#include <armadillo>
#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include <fstream>
using namespace std;
using namespace arma;

Model::Model() {
	// TODO Auto-generated constructor stub

}


Model::Model(string name, double modelCenterX, double modelCenterY, double critRad, double e,  double PA, double mass):
		name(name), modelCenterX(modelCenterX), modelCenterY(modelCenterY), critRad(critRad), e(e), PA(PA), mass(mass) {

}



vector<double> Model::getDeflectionAngle(Conf* conList, int imgX, int imgY) {
	double fX = (imgX -conList->imgXCenter-modelCenterX)*conList->imgRes;   // Model center frame
	double fY = (imgY -conList->imgYCenter-modelCenterY)*conList->imgRes;

	double pfX = (imgX-conList->imgXCenter)*conList->imgRes; 			// image center frame;
	double pfY = (imgY-conList->imgYCenter)*conList->imgRes;

	double pDeltaX = 0;
	double pDeltaY = 0;
	vector<double> srcPos;
	double fDenom, srcX, srcY;

	if(name.compare("PTMASS")==0) {
		fDenom = fX*fX+fY*fY;
		double fMult = critRad*critRad/fDenom;
		pDeltaX = fX*fMult;
		pDeltaY = fY*fMult;

	}

	if(name.compare("SIE")==0) {

		double phi,root1mq,fq,fac,fCore=0,fCosTheta,fSinTheta,x1,y1,deltax1,deltay1;
			//	fX -= g_PixelResn*pLensComp->fParameter[3];
			//	fY -= g_PixelResn*pLensComp->fParameter[4];

				if (fX == 0 && fY == 0)
					pDeltaX = pDeltaY = critRad; // pLensComp->fParameter[0];

				/* pre-calculate constants */
				fCosTheta = cos(PA*M_PI/180);
				fSinTheta = sin(PA*M_PI/180);
				fq = 1-e;
                if (fq>1.0) cout << "Axis ratio should be smaller than 1. " << endl;
                if (fq==1.0) fq = 0.999;

				/* rotate reference frame to x-axis */
				x1 = fX*fCosTheta + fY*fSinTheta;
				y1 = -fX*fSinTheta + fY*fCosTheta;

				root1mq = sqrt(1.0-fq*fq);
				phi = sqrt(fq*fq*(fCore*fCore + x1*x1) + y1*y1);
				/* use sqrt(fq) here not fq. This presevers scale of Einstein Radius */
				fac = critRad*sqrt(fq)/root1mq;
				//cout<< fac << endl;
				deltax1 = fac*atan(root1mq*x1/(phi + fCore));
				deltay1 = fac*lm_arctanh(root1mq*y1/(phi+ fCore*fq*fq));


				pDeltaX = deltax1*fCosTheta - deltay1*fSinTheta;
				pDeltaY = deltay1*fCosTheta + deltax1*fSinTheta;


	}
/*
	if(name.compare("NFW")) {
		 parameters:  mass scale, scale len, ellipticity, orientation_angle
				{
					real_t	fEllip,fCosTheta,fSinTheta,x1,y1,fPhi,fAngRadius,fTempResult,fCosPhi,fSinPhi,fScale;
	                fX -= g_PixelResn*pLensComp->fParameter[4];
	                fY -= g_PixelResn*pLensComp->fParameter[5];
					 pre-calculate constants
					fCosTheta = cos(pLensComp->fParameter[3]*M_PI/180);
					fSinTheta = sin(pLensComp->fParameter[3]*M_PI/180);

					fEllip = pLensComp->fParameter[2];
					fScale = pLensComp->fParameter[1];

	                if (fEllip >= 1.0 || fEllip < 0 || fScale < 0) {
	                    iStatus = LM_BADPROJ;
	                    break;
	                }

					 create elliptical co-ords still in angle units from rotated frame sky coords
					x1 = sqrt(1.0 - fEllip)*(fX*fCosTheta + fY*fSinTheta);
					y1 = sqrt(1.0 + fEllip)*(-fX*fSinTheta + fCosTheta*fY);
					fPhi = atan2(y1,x1);

					 angular radius is in dimensionless units
					fAngRadius = sqrt(x1*x1 + y1*y1)/fScale;

					if (fAngRadius > 0.0) {
						real_t	deflx,defly;

						fCosPhi = cos(fPhi);
						fSinPhi = sin(fPhi);
						 calculate mass. Note that one multiple of fAngRadius has been removed from denominator in next
						 * line so that we don't need to mult again in the following two
						fTempResult = pLensComp->fParameter[0]*fScale*lm_nfw_mass(fAngRadius)/(fAngRadius);
						deflx = sqrt(1.-fEllip)*fTempResult*fCosPhi;
						defly = sqrt(1.+fEllip)*fTempResult*fSinPhi;
						*pDeltaX = (deflx*fCosTheta - defly*fSinTheta);
						*pDeltaY = (deflx*fSinTheta + defly*fCosTheta);
					}
					else {

						iStatus =LM_IGNORE_PROJ;

						*pDeltaX = 0.0;
						*pDeltaY = 0.0;
					}
				}

	}*/


	srcX = (pfX - pDeltaX)/conList->srcRes+conList->srcXCenter;
	srcY = (pfY - pDeltaY)/conList->srcRes+conList->srcYCenter;

	//ofstream debug("debug.txt" , ios::out | ios::app  );

	//debug << imgX << "\t" << imgY  << "\t" <<pDeltaX<< "\t" << pDeltaY << "\t"<<  srcX << "\t" << srcY << endl;
	//debug.close();
	srcPos.push_back(srcX);
	srcPos.push_back(srcY);
	return srcPos;

}


Model::~Model() {
	// TODO Auto-generated destructor stub
}



void Model::updatePosMapping(Image* image, Conf* conList) {

	length = conList->length;
	vector<double> srcPos;
	for(int i=0; i<length; ++i) {
		int imgX = image->xList[i];
		int imgY = image->yList[i];
		srcPos = getDeflectionAngle(conList,imgX, imgY);
		srcPosXList.push_back(srcPos[0]);
		srcPosYList.push_back(srcPos[1]);

		posMap[make_pair(imgX, imgY)] = i;
	}

}

sp_mat Model::buildLensMatrix(Image* dataImage,  Conf* constList) {

	sp_mat L(constList->length,constList->length);
	map<pair<int, int>,int>::iterator left, right, up, down;
	double Ax=0.0, Ay=0.0, Bx=0.0, By=0.0, Cx=0.0, Cy=0.0, Px=0.0, Py=0.0;
	vector<double> w;
	//int counter;
	for (int i=0; i<constList->length; ++i) {
		if(dataImage->type[i]==1)
			L(i,i)=1;
		if(dataImage->type[i]==0) {
			left  = posMap.find(make_pair(dataImage->xList[i]-1, dataImage->yList[i]));
			right = posMap.find(make_pair(dataImage->xList[i]+1, dataImage->yList[i]));
			up    = posMap.find(make_pair(dataImage->xList[i], dataImage->yList[i]+1));
			down  = posMap.find(make_pair(dataImage->xList[i], dataImage->yList[i]-1));

			//counter = ((left!=posMap.end()) + (up!=posMap.end()) + (down!=posMap.end()) + (right!=posMap.end()));
			//cout << counter << endl;

			if(left!=posMap.end() && up!=posMap.end() && down!=posMap.end()) {
				int iLeft = left->second;
				int iUp   = up  ->second;
				int iDown = down->second;

				Ax = srcPosXList[iLeft];
				Ay = srcPosYList[iLeft];
				Bx = srcPosXList[iUp];
				By = srcPosYList[iUp];
				Cx = srcPosXList[iDown];
				Cy = srcPosYList[iDown];
				Px = srcPosXList[i];
				Py = srcPosYList[i];

				w = getTriWeight( Ax,  Ay,  Bx,  By,  Cx,  Cy,  Px,  Py);

				L(i, iLeft) = w[0];
				L(i, iUp)   = w[1];
				L(i, iDown) = w[2];
			}
			else L(i, i) = 1;
		}
	}
	return L;

}
