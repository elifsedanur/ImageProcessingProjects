#pragma once
class Object
{
private:
	int rowStart;
	int rowFinish;
	int colStart;
	int colFinish;
	unsigned char* object;
	double moments[7];
	double centroidX;
	double centroidY;

public:
	Object(int rows,int rowf,int cols,int colf,unsigned char* obj);
	~Object();
	/*double rawMoment(int p, int q);
	void calculateCentroid();
	double centralMoment(int p, int q);
	double normalizedMoment(int p, int q);
	void invariantMoments();
	double* getMoments();
	int getRows();
	int getRowf();
	int getCols();
	int getColf();*/

};

