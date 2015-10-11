#include <iostream>
#include <string>
#include <map>
#include "Image.h"
#include "commons.h"
#include "fitsio.h"
#include "Model.h"
#include <armadillo>
using namespace std;
using namespace arma;


Image* prepare(string imgFileName, string regionFileName, string varFileName) {
	Image *dataImage = new Image(imgFileName);
	Image *varImage = new Image(varFileName);
	dataImage->updateFilterImage(regionFileName);
	dataImage->writeFilterImage("filteredImage.fits");
	sp_mat C= varImage->getVarMatrix(regionFileName);

	return dataImage;
}




int main() {

	string imageFileName="jun_image.fits";
	string regionFileName = "mask_l.reg";
	string confFileName = "conf.txt";
	string varFileName = "jun_var.fits";
	Image* dataImage = prepare( imageFileName, regionFileName, varFileName);
	Conf *conList = new Conf(dataImage, confFileName);

	vec d = dataImage->getMatrixD();

	//cout << d << endl;






	// output source Image:

	vector<double> srcX(conList->length, 0);
	vector<double> srcY(conList->length, 0);

	Model LM_PTMASS("PTMASS", 0, 0, 5.7, 0,0,0 );

	map<pair<int, int>,int> posMap = LM_PTMASS.createPosMapping(dataImage, &srcX, &srcY, conList);

	vector<double> srcBriList(conList->length, 0);
	for (int i=0; i<conList->length; ++i) {
		srcBriList[i] = d[i];
	}


	conList->printConfList();

	Image srcImg(srcX, srcY, &srcBriList, conList->srcSize[0], conList->srcSize[1], conList->bitpix);

	//srcImg.printImageInfo(1,1, 100, 100);
	dataImage->writeToFile("test.fits");
	srcImg.writeToFile("src.fits");



	return 0;


}
