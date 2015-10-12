#include <iostream>
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

	cout << "critRad \t chi2 \t" << endl;
	for (int i=0; i<50; ++i) {
		double critRad = 3.7+i*0.05;
		double e = 0.01 + i*0.01;
		double PA = i*3;

		Model *model = new Model("SIE", 0, 0, 6.16, 0.01, PA,0);
		model->updatePosMapping(dataImage, conList);

		sp_mat L = model->buildLensMatrix(dataImage, conList);

		double chi2 = getPenalty(&L , &d,  &d, &C);

		cout << PA << "\t" << chi2 << endl;

	}

	return 0;


}
