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


int main() {

	/* prepare */
	string confFileName = "conf.txt";
	map<string, string> mapConf = parseConfigure(confFileName);
	Image* dataImage = new Image(mapConf["imageFileName"]);
	dataImage->updateFilterImage(mapConf["regionFileName"]);
	Image* varImage = new Image(mapConf["varFileName"]);
	Conf *conList = new Conf(dataImage, mapConf);
	dataImage->updateGridPointType();
	//dataImage->writeFilterImage("filteredImage.fits");
	sp_mat C = varImage->getVarMatrix(mapConf["regionFileName"]);
	vec d =dataImage->getMatrixD();

	/* Now these are ready:
		dataImag, conList, C, d
	*/

	// Build Matrix L;

	/*vector<double> srcBriList(conList->length, 0);
		for (int i=0; i<conList->length; ++i) {
			srcBriList[i] = d[i];
		}



		conList->printConfList();

		Image* srcImg = new Image(LM_PTMASS.srcPosXList, LM_PTMASS.srcPosYList, &srcBriList, conList->srcSize[0], conList->srcSize[1], conList->bitpix);

		//srcImg.printImageInfo(1,1, 100, 100);
		dataImage->writeToFile("test.fits");
		srcImg->writeToFile("src.fits");

	*/


	ofstream f;
	f.open("list.txt");
	cout << "critRad \t chi2 \t" << endl;
	double minChi2 = 10e6;
	double minE  = 10e6;
	double minCritR = 10e6;
	for (int i=0; i<100; ++i) {
		//for(int j=0; j<40; ++j) {
			double critRad = 3.7+i*0.05;
			double e = 0.2 ;//  + j*0.02;
			double PA = 0 ; //2*j+30;

			Model *model = new Model("SIE", 0, 0, critRad, e, PA,0);
			model->updatePosMapping(dataImage, conList);

			sp_mat L = model->buildLensMatrix(dataImage, conList);

			double chi2 = getPenalty(&L , &d,  &d, &C);
			if(chi2 < minChi2)  {
				minChi2 = chi2;
				minE = e;
				minCritR = critRad;
			}

			f <<e << "\t" <<  critRad << "\t" << chi2 << endl;

		//}
	}
	cout <<minE << "\t" <<  minCritR << "\t" << minChi2 << endl;
	f.close();
	return 0;


}
