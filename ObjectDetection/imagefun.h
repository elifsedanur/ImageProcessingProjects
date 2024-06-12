#pragma once
#include <iostream>
#include <vector>
using namespace std;

void RGBtoGray(unsigned char* input, unsigned char* grayimage, int rows, int cols);
int* Histogram(unsigned char* inp, int h, int w);
unsigned char* kMeans(unsigned char* inp, int r, int c, int* hist, float k1, float k2);
unsigned char* Dilation(unsigned char* inp, int r, int c);
std::vector<int> Labeling(unsigned char* inp, int r, int c, unsigned char* out);
unsigned char* Framing(unsigned char* inp, int r, int c, int tag);
double* Moments1(unsigned char* object, int r, int c);

double calculateDistance(double*, double*);
