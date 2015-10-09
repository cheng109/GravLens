#include <iostream>
#include <string>
#include "commons.h"
#include "Image.h"
#include "fitsio.h"
using namespace std;


int main() {
	
	
	cout<<"hello world" << endl;
	string imgFileNameS = "sample.fits";
	//char* filename = "sample.fits";
	//readimage(filename);
	Image dataImage(imgFileNameS);



	dataImage.printImageInfo(5, 5);

	string outFileName = "output.fits";
	dataImage.writeToFile(outFileName);

	string regionFileName = "mask.reg";


	//parseReagionFile( regionFileName);

}
