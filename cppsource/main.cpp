#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "Image.h"
#include "commons.h"
#include "fitsio.h"
#include "Model.h"
#include <armadillo>
#include <tuple>
#include <map>
using namespace std;
using namespace arma;
//#include "gnuplot-iostream.h"

int main() {

	// prepare
	string confFileName = "conf.txt";
	map<string, string> mapConf = parseConfigure(confFileName);
	Image* dataImage = new Image(mapConf["imageFileName"]);
	dataImage->updateFilterImage(mapConf["regionFileName"]);
	dataImage->updateVarList(2.2, 0.1);
	sp_mat invC = dataImage->getVarMatrix();

	Conf *conList = new Conf(dataImage, mapConf);

	dataImage->updateGridPointType();

	dataImage->writeFilterImage("filteredImage.fits");

	vec d =dataImage->getMatrixD();
	// Now these are ready:
	//	dataImag, conList, C, d

	vector<double> srcBriList(conList->length, 0);
	for (int i=0; i<conList->length; ++i) {
		srcBriList[i] = d[i];
	}

	ParaList para("PTMASS", 0, 0, 5.69, 0.0, 0.0, 0.0, 0.0);


	Model *model = new Model(conList, para);
	model->updatePosMapping(dataImage, conList);
	Image* srcImg = new Image(model->srcPosXList, model->srcPosYList, &srcBriList, conList->srcSize[0], conList->srcSize[1], conList->bitpix);
	//srcImg.printImageInfo(1,1, 100, 100);
	dataImage->writeToFile("test.fits");
	srcImg->writeToFile("src.fits");


	model->Logging(dataImage, conList, "log.txt");

/*
	double minChi2 = 10e6;
	double minE  = 10e6;
	double minCritR = 10e6;
	for (int i=0; i<1; ++i) {
		//for(int j=0; j<40; ++j) {
			double critRad = 3.95+i*0.05;
			double e = 0.0 ;//  + j*0.02;
			double PA = 0 ; //2*j+30;
			//cout << i << endl;
			Model *model = new Model(conList, "PTMASS", 0, 0, critRad, e, PA,0);
			model->updatePosMapping(dataImage, conList);
			sp_mat L = model->buildLensMatrix(dataImage, conList);
			double chi2 = getPenalty(&L , &d,  &d, &invC);
			if(chi2 < minChi2)  {
				minChi2 = chi2;
				minE = e;
				minCritR = critRad;
			}
			//y.push_back(chi2);
			cout <<e << "\t" <<  critRad << "\t" << chi2 << endl;
		//}
	}
	cout <<minE << "\t" <<  minCritR << "\t" << minChi2 << endl;
	//gp.send(y);
	conList->printConfList();*/
	return 0;
}

