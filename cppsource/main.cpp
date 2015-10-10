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


Image* prepare(string imgFileName, string regionFileName) {
	Image *dataImage = new Image(imgFileName);
	dataImage->updateFilterImage(regionFileName);
	dataImage->writeFilterImage("filteredImage.fits");


	return dataImage;

}





int main() {

	string imageFileName="jun_image.fits";
	string regionFileName = "mask.reg";


	Image* dataImage = prepare( imageFileName, regionFileName);
	Const *conList = new Const(dataImage);

	vec d = dataImage->getMatrixD();

	//cout << d << endl;






	// output source Image:

	vector<double> srcX(conList->length, 0);
	vector<double> srcY(conList->length, 0);

	Model LM_PTMASS("PTMASS", 0, 0, 5.7, 0,0,0 );

	map<pair<int, int>,int> posMap = LM_PTMASS.createPosMapping(dataImage, &srcX, &srcY, conList);
	for(int i=0; i<conList->length; ++i) {
	if(srcX[i]>60) cout << srcX[i] << "\t" << srcY[i] << endl; }
	vector<double> srcBriList(conList->length, 0);
	for (int i=0; i<conList->length; ++i) {
		srcBriList[i] = d[i];
	}

	Image srcImg(srcX, srcY, &srcBriList, conList->srcSize[0], conList->srcSize[1], conList->bitpix);
	//srcImg.printImageInfo(1,1, 100, 100);
	dataImage->writeToFile("test.fits");
	srcImg.writeToFile("src.fits");

	conList->printConstList();

	return 0;


}
