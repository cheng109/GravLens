#include <iostream>
#include <string>
#include "commons.h"
#include "Image.h"
#include "fitsio.h"
using namespace std;



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
	Const *constList = new Const(dataImage);
	constList->printConstList();





	//dataImage.printImageInfo(22,14,25,21);


	//dataImage.writeToFile("output.fits");



	vector<double> xpos;
	vector<double> ypos;


	size_t n = parseReagionFile("mask.reg", &xpos, &ypos);







}
