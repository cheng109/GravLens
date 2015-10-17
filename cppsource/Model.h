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
	vector<double> something;
	sp_mat L;
	sp_mat M;
	sp_mat r;
	sp_mat s0;


	sp_mat Ds;
	sp_mat Dphi;
	sp_mat Hs1;
	sp_mat Hs2;
	sp_mat Hphi;
	sp_mat HtH;
	sp_mat HphiH;
	double lambdaS;
	double lambdaPhi;
	sp_mat RtR;

	double penalty;
	vector<vector<normVec> > normV;
	vector<normVec> meanNormV;
	map<pair<int, int>,int> posMap;



public:
	Model();
	Model(Conf* conList, ParaList param);
	vector<double> getDeflectionAngle(Conf* conList, int imgX, int imgY);
	void updatePosMapping(Image* image,  Conf* conList);
	sp_mat buildLensMatrix(Image* dataImage,  Conf* constList);
	void updateGradient(Image* dataImage);
	void Logging(Image* dataImage, Conf* conList, string outFileName);
	void updateRegularMatrix();
	double getPenalty();
	virtual ~Model();
};

#endif /* MODEL_H_ */
