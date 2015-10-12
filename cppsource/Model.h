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
	double critRad;
	double e;
	double q;
	double PA;
	double mass;

	int length;
public:
	vector<double> srcPosXList;
	vector<double> srcPosYList;

	map<pair<int, int>,int> posMap;


public:
	Model();
	Model(string name,double modelCenterX, double modelCenterY, double critR, double e, double PA, double mass);
	vector<double> getDeflectionAngle(Conf* conList, int imgX, int imgY);
	void updatePosMapping(Image* image,  Conf* conList);
	sp_mat buildLensMatrix(Image* dataImage,  Conf* constList);
	virtual ~Model();
};

#endif /* MODEL_H_ */
