/*
 * Model.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: juncheng
 */

#include "Model.h"

Model::Model() {
	// TODO Auto-generated constructor stub

}



void Model::getDeflectionAnle(Const* conList, size_t imgX, size_t imgY, double *dx, double *dy) {
	if(name.compare("PTMASS")==0) {
		double fX = (imgX -conList->imgXCenter-modelCenterX)*conList->imgRes;
		double fY = (imgY -conList->imgYCenter-modelCenterY)*conList->imgRes;
	}


}



Model::~Model() {
	// TODO Auto-generated destructor stub
}

