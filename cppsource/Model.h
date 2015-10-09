/*
 * Model.h
 *
 *  Created on: Oct 9, 2015
 *      Author: juncheng
 */

#ifndef MODEL_H_
#define MODEL_H_
#include <string>
#include <vector>
#include "commons.h"
using namespace std;

class Model {

	string name;
	double modelCenterX;
	double modelCenterY;


public:
	Model();
	void getDeflectionAnle(Const* conList, size_t imgX, size_t imgY, double *dx, double *dy);
	virtual ~Model();
};

#endif /* MODEL_H_ */
