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


Model::Model(Conf* conList, ParaList param):
		param(param), length(conList->length),
		L(length,length),
		Ds(length, 2*length),
		Dphi(2*length,length),
		Hs1(length, length),
		Hs2(length, length),
		Hphi(length, length),
		HtH(length, length),
		HphiH(length, length),
		RtR(2*length, 2*length) {
	// initial s;

	normVec n(0, 0, 0);
	vector<normVec> temp;
	temp.push_back(n);
	//sp_mat L(constList->length,constList->length);
	//normV.resize(conList->length);

	for(int i=0; i<conList->length; ++i) {
		s.push_back(0);
		normV.push_back(temp);
		dSy1.push_back(0);
		dSy2.push_back(0);
	}
}



vector<double> Model::getDeflectionAngle(Conf* conList, int imgX, int imgY) {
	double fX = (imgX -conList->imgXCenter-param.centerX)*conList->imgRes;   // Model center frame
	double fY = (imgY -conList->imgYCenter-param.centerY)*conList->imgRes;

	double pfX = (imgX-conList->imgXCenter)*conList->imgRes; 			// image center frame;
	double pfY = (imgY-conList->imgYCenter)*conList->imgRes;

	double pDeltaX = 0;
	double pDeltaY = 0;
	vector<double> srcPos;
	double fDenom, srcX, srcY;



	if(param.name.compare("PTMASS")==0) {
		fDenom = fX*fX+fY*fY;
		double fMult = param.critRad*param.critRad/fDenom;
		pDeltaX = fX*fMult;
		pDeltaY = fY*fMult;

	}

	/*if(name.compare("SIE")==0) {

		double phi,root1mq,fq,fac,fCore=0,fCosTheta,fSinTheta,x1,y1,deltax1,deltay1;
			//	fX -= g_PixelResn*pLensComp->fParameter[3];
			//	fY -= g_PixelResn*pLensComp->fParameter[4];

				if (fX == 0 && fY == 0)
					pDeltaX = pDeltaY = critRad; // pLensComp->fParameter[0];

				 pre-calculate constants
				fCosTheta = cos(PA*M_PI/180);
				fSinTheta = sin(PA*M_PI/180);
				fq = 1-e;
                if (fq>1.0) cout << "Axis ratio should be smaller than 1. " << endl;
                if (fq==1.0) fq = 0.999;

				 rotate reference frame to x-axis
				x1 = fX*fCosTheta + fY*fSinTheta;
				y1 = -fX*fSinTheta + fY*fCosTheta;

				root1mq = sqrt(1.0-fq*fq);
				phi = sqrt(fq*fq*(fCore*fCore + x1*x1) + y1*y1);
				 use sqrt(fq) here not fq. This presevers scale of Einstein Radius
				fac = critRad*sqrt(fq)/root1mq;
				//cout<< fac << endl;
				deltax1 = fac*atan(root1mq*x1/(phi + fCore));
				deltay1 = fac*lm_arctanh(root1mq*y1/(phi+ fCore*fq*fq));


				pDeltaX = deltax1*fCosTheta - deltay1*fSinTheta;
				pDeltaY = deltay1*fCosTheta + deltax1*fSinTheta;


	}*/
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


	map<pair<int, int>,int>::iterator left, right, up, down;
	vector<double> w;
	//int counter;
	for (int i=0; i<constList->length; ++i) {
		left  = posMap.find(make_pair(dataImage->xList[i]-1, dataImage->yList[i]));
		right = posMap.find(make_pair(dataImage->xList[i]+1, dataImage->yList[i]));
		up    = posMap.find(make_pair(dataImage->xList[i], dataImage->yList[i]+1));
		down  = posMap.find(make_pair(dataImage->xList[i], dataImage->yList[i]-1));

		if(left!=posMap.end() && up!=posMap.end() && down!=posMap.end() && right!=posMap.end()) {
			int iLeft = left->second;
			int iUp   = up  ->second;
			int iDown = down->second;
			int iRight= right->second;
			Point A(srcPosXList[iLeft], srcPosYList[iLeft], s[iLeft]);
			Point B(srcPosXList[iUp], srcPosYList[iUp], s[iUp]);
			Point C(srcPosXList[i], srcPosYList[i], s[i]);
			Point D(srcPosXList[iDown], srcPosYList[iDown], s[iDown]);
			Point E(srcPosXList[iRight], srcPosYList[iRight], s[iRight]);
			vector<double> w = getPentWeigth(A, B, C, D, E);
			Hs1(i, iLeft) 	= w[0];
			Hs1(i, iUp) 	= w[1];
			Hs1(i, i) 		= w[2];
			Hs1(i, iDown) 	= w[3];
			Hs1(i, iRight) 	= w[4];

			Hs2(i, iLeft) 	= w[5];
			Hs2(i, iUp) 	= w[6];
			Hs2(i, i) 		= w[7];
			Hs2(i, iDown) 	= w[8];
			Hs2(i, iRight) 	= w[9];
		}
		if(dataImage->type[i]==1) {
			L(i,i)=1;

		}
		if(dataImage->type[i]==0) {
			if(left!=posMap.end() && up!=posMap.end() && down!=posMap.end()) {
				int iLeft = left->second;
				int iUp   = up  ->second;
				int iDown = down->second;

				Point A(srcPosXList[iLeft], srcPosYList[iLeft], s[iLeft]);
				Point B(srcPosXList[iUp], srcPosYList[iUp], s[iUp]);
				Point C(srcPosXList[iDown], srcPosYList[iDown], s[iDown]);
				Point P(srcPosXList[i], srcPosYList[i], s[i]);

				w = getTriWeight( A, B, C, P);
				L(i, iLeft) = w[0];
				L(i, iUp)   = w[1];
				L(i, iDown) = w[2];
				normVec n = getNormVector(A, B, C);
				normV[iLeft].push_back(n);
				normV[iUp].push_back(n);
				normV[iDown].push_back(n);

			}
			else L(i, i) = 1;
		}
	}
	HtH = Hs1.t()*Hs1 + Hs2.t()*Hs2;
	for(sp_mat::iterator iter=HtH.begin(); iter!=HtH.end(); ++iter) {

		RtR(iter.row(), iter.col()) = *iter;
	}
	mat RR(RtR);
	cout << RR << endl;
	return L;
}

void Model::updateGradient(Image* dataImage) {
	vec vSy1(length);  vSy1.fill(0);
	vec vSy2(length);  vSy2.fill(0);
	for (int i=0; i<length; ++i) {
		if(dataImage->type[i]==1) {
			normVec mean =  meanNormVector(normV[i]);
			if (mean.n2==0) {
				vSy1[i] = 0;
				vSy2[i] = 0;
			}
			else {
				vSy1[i] = -mean.n0/mean.n2;
				vSy2[i] = -mean.n1/mean.n2;
			}
		}
	}
	vSy1 = L*vSy1;
	vSy2 = L*vSy2;
	//update Ds and Dphi
	for(int i=0; i<length; ++i){
		for(int j=0; j<2*length; ++j) {
			Ds(i, 2*i+0) = vSy1(i);
			Ds(i, 2*i+1) = vSy2(i);
		}
	}

	sp_mat extension = -1*Ds*Dphi;
	mat phi(length, 1);
	M = join_rows(L, extension);
	mat sMatrix(length, 1);
	for(int i=0; i<length; ++i) {
		sMatrix(i,0) = s[i];

	}
	r = join_cols(sMatrix,  phi);

}

double Model::getPenalty() {


	return 0;
}

void Model::Logging(Image* dataImage, Conf* conList, string outFileName) {
	ofstream f(outFileName);
	string tab = "\t";
	string entry;
	f << "#1 index\n" << "#2 imgX\n" << "#3 imgY\n" << "#4 imgBri\n";
	f << "#5 srcX\n"  << "#6 srcY\n" << "#7 mapIndex" << "#8 ";
	for(int i=0; i<conList->length; ++i) {
		entry =  to_string(dataImage->iList[i]) + "\t"
				+to_string(dataImage->xList[i]) + "\t"
				+to_string(dataImage->yList[i]) + "\t"
				+to_string(dataImage->dataList[i]) + "\t"
				+to_string(srcPosXList[i]) + "\t"
				+to_string(srcPosYList[i]) + "\t"
				+to_string(posMap[make_pair(dataImage->xList[i], dataImage->yList[i])]) + "\t" ;  //  Matched position
				//+to_string()

		f << entry << endl ;
	}

	f.close();


}

