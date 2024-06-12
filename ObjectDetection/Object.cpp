#include "Object.h"
#include <math.h>
#include <iostream>
using namespace std;

Object::Object(int rows, int rowf, int cols, int colf, unsigned char* obj) {
	this->rowStart = rows;
	this->rowFinish = rowf;
	this->colStart = cols;
	this->colFinish = colf;
	this->object = new unsigned char[(rowFinish - rowStart + 1) * (colFinish - colStart + 1)];

	/*calculateCentroid();
	invariantMoments();*/
}
Object::~Object() {
	delete[] object;
}
//double Object::rawMoment(int p, int q) {
//	double sum = 0.0;
//	int row = (rowFinish - rowStart + 1);
//	int col = (colFinish - colStart + 1);
//	for (int i = 0; i < row; i++) {
//		for (int j = 0; j < col; j++) {
//			sum += pow(i, p) * pow(j, q) * object[i * col + j];
//		}
//	}
//	/*for (int i = 0; i < row; i++) {
//		for (int j = 0; j < col; j++) {
//			rawMoment += image[i * col + j] * pow(i, p) * pow(j, q);
//		}
//	}*/
//	cout << "raw moment(" << p << "," << q << "): " << sum << endl;
//	return sum;
//}
//void Object::calculateCentroid() {
//	this->centroidX = rawMoment(1, 0) / rawMoment(0, 0);
//	cout << "r: " << centroidX << endl;
//	this->centroidY = rawMoment(0, 1) / rawMoment(0, 0);
//	cout << "c: " << centroidY << endl;
//}
//double Object::centralMoment(int p, int q) {
//	double sum = 0.0;
//	int row = (rowFinish - rowStart + 1);
//	int col = (colFinish - colStart + 1);
//	for (int i = 0; i < row; i++) {
//		for (int j = 0; j < col; j++) {
//			sum += pow(i-centroidX, p) * pow(j-centroidY, q) * object[i * col + j];
//		}
//	}
//	return sum;
//}
//double Object::normalizedMoment(int p, int q) {
//	int row = (rowFinish - rowStart + 1);
//	int col = (colFinish - colStart + 1);
//	int y = (p + q) / 2 + 1;
//
//	double deg = centralMoment(p, q) / pow(centralMoment(0, 0), y);
//	return deg;
//}
//
//void Object::invariantMoments() {
//	moments[0] = normalizedMoment(2, 0) + normalizedMoment(0, 2);
//	moments[1] = pow(normalizedMoment(2, 0) - normalizedMoment(0, 2), 2) + 4 * pow(normalizedMoment(1, 1), 2);
//	moments[2] = pow(normalizedMoment(3, 0) - 3*normalizedMoment(1, 2), 2) + 
//		pow(3 * normalizedMoment(2, 1) - normalizedMoment(0, 3), 2);
//	moments[3] = pow(normalizedMoment(3, 0) + normalizedMoment(1, 2), 2) + pow(normalizedMoment(2, 1) + normalizedMoment(0, 3), 2);
//	moments[4] = (normalizedMoment(3, 0) - 3 * normalizedMoment(1, 2)) * (normalizedMoment(3, 0) + normalizedMoment(1, 2)) *
//		(pow(normalizedMoment(3, 0) + normalizedMoment(1, 2), 2) - 3 * pow(normalizedMoment(2, 1) + normalizedMoment(0, 3), 2)) +
//		(3 * normalizedMoment(2, 1) - normalizedMoment(0, 3)) * (normalizedMoment(2, 1) + normalizedMoment(0, 3)) *
//		(3 * pow(normalizedMoment(3, 0) + normalizedMoment(1, 2), 2) - pow(normalizedMoment(2, 1) + normalizedMoment(0, 3), 2));
//	moments[5] = (normalizedMoment(2, 0) - normalizedMoment(0, 2)) * (pow(normalizedMoment(3, 0) + normalizedMoment(1, 2), 2) - pow(normalizedMoment(2, 1) + normalizedMoment(0, 3), 2)) +
//		4 * normalizedMoment(1, 1) * (normalizedMoment(2, 1) + normalizedMoment(0, 3))*(normalizedMoment(3,0)+normalizedMoment(1,2));
//	moments[6] = (3*normalizedMoment(2, 1) - normalizedMoment(0, 3)) * (normalizedMoment(3, 0) + normalizedMoment(1, 2)) *
//		(pow(normalizedMoment(3, 0) + normalizedMoment(1, 2), 2) - 3 * pow(normalizedMoment(2, 1) + normalizedMoment(0, 3), 2)) -
//		(3 * normalizedMoment(1,2) - normalizedMoment(3,0)) * (normalizedMoment(2, 1) + normalizedMoment(0, 3)) *
//		(3 * pow(normalizedMoment(3, 0) + normalizedMoment(1, 2), 2) - pow(normalizedMoment(2, 1) + normalizedMoment(0, 3), 2));
//
//}
//
//double* Object::getMoments() {
//	
//	return this->moments;
//}
//int Object::getRows() { return this->rowStart; }
//int Object::getRowf() { return this->rowFinish; }
//int Object::getCols() { return this->colStart; }
//int Object::getColf() { return this->colFinish; }




