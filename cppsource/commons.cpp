/*
 * commons.cpp
 *
 *  Created on: Sep 30, 2015
 *      Author: cheng109
 */
#include "commons.h"
#include "fitsio.h"
#include <string>
#include "Image.h"
#include <iostream>
#include <fstream>
#include <sstream>

#define buffsize 1000

using namespace std;


void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/
    if (status)
    {
       fits_report_error(stderr, status); /* print error report */
       exit( status );    /* terminate the program, returning error status */
    }
    return;
}

bool pnpoly(size_t nvert, vector<double> *vertx, vector<double> *verty, double testx, double testy)
{
  int i, j;
  bool c = FALSE;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty->at(i)>testy) != (verty->at(j)>testy)) &&
     (testx < (vertx->at(j)-vertx->at(i)) * (testy-verty->at(i)) / (verty->at(j)-verty->at(i)) + vertx->at(i)) )
       c = !c;
  }
  return c;
}


size_t parseReagionFile(string regionFileName, vector<double> *xpos, vector<double> *ypos) {

	ifstream regionFile(regionFileName.c_str());
	string line, token;
	size_t pos=0;
	//vector<double> xpos;
	//vector<double> ypos;

	while (getline(regionFile, line)) {
		if (line[0]!='#' && line.substr(0, 6)!="global" && line.substr(0, 5)!="image") {
			size_t pos1 = line.find_first_of("(", pos);
			size_t pos2 = line.find_first_of(")", pos);

			istringstream ss(line.substr(pos1+1, pos2-pos1));
			int flag=-1;
			while(getline(ss, token, ',' )){
				if(flag<0) xpos->push_back(stod(token));
				if(flag>0) ypos->push_back(stod(token));
				cout << token << endl;
				flag = (-1)*flag;
			};
		}
		if(xpos->size()!=ypos->size())
			cout << "Error when reading region File!" << endl;



	}
	return xpos->size();

}







