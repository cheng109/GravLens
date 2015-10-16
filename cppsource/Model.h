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

struct ParaList{
	string name;
	double centerX;
	double centerY;
	double critRad;
	double e;
	double q;
	double PA;
	double mass;
	double coreRad;
	ParaList() {};
	ParaList(string name, double x, double y, double critRad, double e, double PA, double mass, double coreRad):
		name(name), centerX(x), centerY(y), critRad(critRad), e(e),  PA(PA), mass(mass), coreRad(coreRad) {};

};


class Model {
	int length;
	ParaList param;

public:
	vector<double> srcPosXList;
	vector<double> srcPosYList;
	vector<double> dSy1; 
	vector<double> dSy2; 
	vector<double> s; 

	vector<vector<normVec> > normV;
	map<pair<int, int>,int> posMap;


public:
	Model();
	Model(Conf* conList, ParaList param);
	vector<double> getDeflectionAngle(Conf* conList, int imgX, int imgY);
	void updatePosMapping(Image* image,  Conf* conList);
	sp_mat buildLensMatrix(Image* dataImage,  Conf* constList);
	void updateNorm(Image* dataImage);
	void Logging(Image* dataImage, Conf* conList, string outFileName);
	virtual ~Model();
};

#endif /* MODEL_H_ */
