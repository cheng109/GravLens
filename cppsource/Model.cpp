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


void Model::getDeflectionAnle(Const* const, size_t imgX, size_t imgY, double *dx, double *dy) {
	if(name.compare("PTMASS")==0) {
		double fX = (imgX -const->imgXCenter-modelCenterX)*const->imgRes;
		double fY = (imgY -const->imgYCenter-modelCenterY)*const->imgRes;


	}

}


Model::~Model() {
	// TODO Auto-generated destructor stub
}

