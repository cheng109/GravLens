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
#include "Image.h"
#include <map>
using namespace std;

class Model {

	string name;
	double modelCenterX;
	double modelCenterY;
	double critR;
	double e;
	double q;
	double PA;
	double mass;


public:
	Model();
	Model(string name,double modelCenterX, double modelCenterY, double critR, double e, double PA, double mass);
	void getDeflectionAngle(Conf* conList, int imgX, int imgY, double *srcX, double *srcY);
	map<pair<int, int>,int> createPosMapping(Image* image, vector<double>* srcX,vector<double>* srcY,  Conf* conList);
	virtual ~Model();
};

#endif /* MODEL_H_ */
