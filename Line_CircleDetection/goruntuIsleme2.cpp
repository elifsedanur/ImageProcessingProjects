// goruntuIsleme2.cpp : Bu dosya 'main' işlevi içeriyor. Program yürütme orada başlayıp biter.
//
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include <stdlib.h>
#include <string.h>
#include <cstdint>
#include <iomanip>
//#include "imagefun.cpp"


#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>

using namespace std;
using namespace cv;
constexpr double PI = 3.14159265358979323846;
int* gradientDirection;
int rows;
int cols;
#define HORIZONTAL 0
#define VERTICAL 1
#define POSITIVE_DIAGONAL 2
#define NEGATIVE_DIAGONAL 3
Mat image;
struct CircleCenter {
	int x;
	int y;
};


void RGBtoGray(unsigned char* input, unsigned char* grayimage, int rows, int cols) {

	for (int r = 0; r < rows; r++) {
		for (int c = 0; c < cols; c++) {
			grayimage[r * cols + c] = input[r * cols * 3 + c * 3 + 2] * 0.299 + input[r * cols * 3 + c * 3 + 1] * 0.587 + input[r * cols * 3 + c * 3] * 0.114;
		}
	}
}
int* GradientImage(unsigned char* input,int &rows,int &cols) {
	unsigned char* ninput = new unsigned char[(rows + 2) * (cols + 2)];
	int rinp = rows;
	int cinp = cols;
	int sobelX[3][3] = { {-1, 0, 1},
						{-2, 0, 2},
						{-1, 0, 1} };

	int sobelY[3][3] = { {-1, -2, -1},
						 { 0,  0,  0},
						 { 1,  2,  1} };


	int orow = (rows - 3) / 1 + 1;
	int ocol = (cols - 3) / 1 + 1;
	//result 0-1020 arasında değer alabilir.
	int* result = new int[orow * ocol];
	gradientDirection = new int[orow * ocol];
	double angle;
	int temp;
	for (int r = 0; r < orow; r++) {
		for (int c = 0; c < ocol; c++) {
			int gradientX = 0;
			int gradientY = 0;
			for (int m = 0; m < 3; m++) {
				for (int n = 0; n < 3; n++) {
					gradientX += input[(r + m) * cols + c + n] * sobelX[m][n];
					gradientY += input[(r + m) * cols + c + n] * sobelY[m][n];
					//gradientX += ninput[(r + m) * (cinp + 2) + c + n] * ker[m][n];
				}
			}
			
			result[r * ocol + c] = abs(gradientX) + abs(gradientY);

			//atan2 radyan cinsinden değer döndürür 1 radyan 180/pi
			int theta = (int)atan2(gradientY, gradientX);
			theta = theta * (180.0 / M_PI);
			//cout << "theta: " << theta << endl; //Theta negatif değerlerde alıyor.
			theta = (theta + 360) % 360; //0-360 arasına sınırladık
			gradientDirection[r * ocol + c] = theta;

		}
	}

	rows = orow;
	cols = ocol;
	return result;
}
unsigned char* NormalizeValue(int* image, int width, int height) {
	unsigned char* normalizeImage = new unsigned char[height * width];
	int min = INT_MAX;
	int max = INT_MIN;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (image[i * width + j] < min) {
				min = image[i * width + j];
			}
			if (image[i * width + j] > max) {
				max = image[i * width + j];
			}
		}
	}
	for (int i = 0; i < height * width; ++i) {
		normalizeImage[i] = (unsigned char)((image[i] - min) * 255.0 / (max - min));
	}
	return normalizeImage;
}
int findGradientDirection(int gradientDirection) {
	int theta = gradientDirection;
	//Yönler arasında 45 derece bıraktık
	if ((theta <= 180 / 8 && theta >= 15 * 180 / 8) || (theta >= 7 * 180 / 8 && theta <= 9 * 180 / 8) || theta == 0) {
		return 1;  // "-"
	}
	else if ((theta <= 5 * 180 / 8 && theta >= 3 * 180 / 8) || (theta >= 11 * 180 / 8 && theta <= 13 * 180 / 8)) {
		return 2;  // "|"
	}
	else if ((theta <= 3 * 180 / 8 && theta >= 180 / 8) || (theta >= 9 * 180 / 8 && theta <= 11 * 180 / 8)) {
		return 3;  // "/"
	}
	else if ((theta <= 7 * 180 / 8 && theta >= 5 * 180 / 8) || (theta >= 13 * 180 / 8 && theta <= 15 * 180 / 8)) {
		return 4; // "\"
	}
	else {
		std::cout << "error " << theta << std::endl;
	}

}
void applyNonMaximumSuppression(int* gradientImage, int* outputImage,int rows,int cols) {
	/*int rows = 6;
	int cols = 6;*/
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			int direction = findGradientDirection(gradientDirection[i * cols + j]);

			switch (direction) {
			case 1:

				if (j == 0 && gradientImage[i * cols + j] > gradientImage[i * cols + j + 1]) {
					outputImage[i * cols + j] = gradientImage[i * cols + j];
				}
				else if (j == cols - 1 && gradientImage[i * cols + j] > gradientImage[i * cols + j - 1]) {
					outputImage[i * cols + j] = gradientImage[i * cols + j];
				}
				else if (gradientImage[i * cols + j] > gradientImage[i * cols + j - 1] && gradientImage[i * cols + j] > gradientImage[i * cols + j + 1]) {
					outputImage[i * cols + j] = gradientImage[i * cols + j];
				}
				else {
					outputImage[i * cols + j] = 0;
				}
				break;
			case 2:
				if (i == 0 && gradientImage[i * cols + j] > gradientImage[(i + 1) * cols + j]) {
					outputImage[i * cols + j] = gradientImage[i * cols + j];
				}
				else if (i == rows - 1 && gradientImage[i * cols + j] > gradientImage[(i - 1) * cols + j]) {
					outputImage[i * cols + j] = gradientImage[i * cols + j];
				}
				else if (gradientImage[i * cols + j] > gradientImage[(i - 1) * cols + j] && gradientImage[i * cols + j] > gradientImage[(i + 1) * cols + j]) {
					outputImage[i * cols + j] = gradientImage[i * cols + j];
				}
				else {
					outputImage[i * cols + j] = 0;
				}
				break;
			case 4: //"\"
				if ((i == 0 && j == cols - 1) || (i == rows - 1 && j == 0)) {
					outputImage[i * cols + j] = 0;
				}
				else if (i == 0 || j == 0) {
					if (gradientImage[i * cols + j] > gradientImage[(i + 1) * cols + j + 1]) {
						outputImage[i * cols + j] = gradientImage[i * cols + j];
					}
					else {
						outputImage[i * cols + j] = 0;
					}
				}
				else if (i == rows - 1 || j == cols - 1) {
					if (gradientImage[i * cols + j] > gradientImage[(i - 1) * cols + j - 1]) {
						outputImage[i * cols + j] = gradientImage[i * cols + j];
					}
					else {
						outputImage[i * cols + j] = 0;
					}
				}
				else if (gradientImage[i * cols + j] > gradientImage[(i - 1) * cols + j - 1] && gradientImage[i * cols + j] > gradientImage[(i + 1) * cols + j + 1]) {
					outputImage[i * cols + j] = gradientImage[i * cols + j];
				}
				else {
					outputImage[i * cols + j] = 0;
				}
				break;
			case 3: //"/"
				if ((i == 0 && j == 0) || (i == rows - 1 && j == cols - 1)) {
					outputImage[i * cols + j] = 0;
				}
				else if (i == 0 || j == cols - 1) {
					if (gradientImage[i * cols + j] > gradientImage[(i + 1) * cols + j - 1]) {
						outputImage[i * cols + j] = gradientImage[i * cols + j];
					}
					else {
						outputImage[i * cols + j] = 0;
					}
				}
				else if (i == rows - 1 || j == 0) {
					if (gradientImage[i * cols + j] > gradientImage[(i - 1) * cols + j + 1]) {
						outputImage[i * cols + j] = gradientImage[i * cols + j];
					}
					else {
						outputImage[i * cols + j] = 0;
					}
				}
				else if (gradientImage[i * cols + j] > gradientImage[(i - 1) * cols + j + 1] && gradientImage[i * cols + j] > gradientImage[(i + 1) * cols + j - 1]) {
					outputImage[i * cols + j] = gradientImage[i * cols + j];
				}
				else {
					outputImage[i * cols + j] = 0;
				}
				break;
			default:
				outputImage[i * cols + j] = 0;
				break;
			}
		}
	}
}
int* Histogram(int* inp) {
	int* hist = new int[1021];
	for (int i = 0; i < 1021; i++) {
		hist[i] = 0;
	}
	for (int i = 0; i < rows * cols; i++) {
		hist[inp[i]]++;
	}
	//delete[] hist;
	return hist;
}
void hysteresisThresholding(int* inputImage, unsigned char* outputImage,int rows,int cols) {

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			int index = i * cols + j;
			int pixelValue = inputImage[index];

			if (pixelValue >= 120) {
				//if (pixelValue >= 10) {
				outputImage[index] = 255;
			}
			else if (pixelValue > 90) { //arada kalan değerler için hata yönetimini kolaylaştırmak için bilindik bir değer atandı
				outputImage[index] = 100;
			}
			else {
				outputImage[index] = 0;
			}
		}
	}}
void edgesWithHysteresis(unsigned char* input, unsigned char* output,int rows,int cols) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			/*int i = 0;
			int j = 563;*/
			int pix = i * cols + j;
			if (input[pix] == 255) {
				output[pix] = 255;
			}
			//edge direc. gradient direc'a dik 
			else if (input[pix] == 100) {
				int direction = findGradientDirection(gradientDirection[pix]);
				switch (direction) {
				case 1:
					if (i == 0) {
						if(input[(i + 1) * cols + j] == 255)
							output[i * cols + j] = 255;
						else
							output[i * cols + j] = 0;
					}
					else if (i == rows - 1) {
						if(input[(i - 1) * cols + j] == 255)
							output[i * cols + j] = 255;
						else
							output[i * cols + j] = 0;
					}
					else if (input[(i - 1) * cols + j] == 255 || input[(i + 1) * cols + j] == 255) {
						output[i * cols + j] = 255;
					}
					else {
						output[i * cols + j] = 0;
					}
					break;
				case 2:
					if (j == 0) {
						if(input[i * cols + j + 1] == 255)
							output[i * cols + j] = 255;
						else
							output[i * cols + j] = 0;
					}
					else if (j == cols - 1) {
						if(input[i * cols + j - 1] == 255)
							output[i * cols + j] = 255;
						else output[i * cols + j] = 0;
					}
					else if (input[i * cols + j - 1] == 255 || input[i * cols + j + 1] == 255) {
						output[i * cols + j] = 255;
					}
					else {
						output[i * cols + j] = 0;
					}
					break;
				case 4:
					if ((i == 0 && j == 0) || (i == rows - 1 && j == cols - 1)) {
						output[i * cols + j] = 0;
					}
					else if (i == 0 || j == cols - 1) {
						if (input[(i + 1) * cols + j - 1] == 255) {
							output[i * cols + j] = 255;
						}
						else {
							output[i * cols + j] = 0;
						}
					}
					else if (i == rows - 1 || j == 0) {
						if (input[(i - 1) * cols + j + 1] == 255) {
							output[i * cols + j] = 255;
						}
						else {
							output[i * cols + j] = 0;
						}
					}
					else if (input[(i - 1) * cols + j + 1] == 255 || input[(i + 1) * cols + j - 1] == 255) {
						output[i * cols + j] = 255;
					}
					else {
						output[i * cols + j] = 0;
					}
					break;
				case 3:
					if ((i == 0 && j == cols - 1) || (i == rows - 1 && j == 0)) {
						output[i * cols + j] = 0;
					}
					else if (i == 0 || j == 0) {
						if (input[(i + 1) * cols + j + 1] == 255) {
							output[i * cols + j] = 255;
						}
						else {
							output[i * cols + j] = 0;
						}
					}
					else if (i == rows - 1 || j == cols - 1) {
						if (input[(i - 1) * cols + j - 1] == 255) {
							output[i * cols + j] = 255;
						}
						else {
							output[i * cols + j] = 0;
						}
					}
					else if (input[(i - 1) * cols + j - 1] == 255 || input[(i + 1) * cols + j + 1] == 255) {
						output[i * cols + j] = 255;
					}
					else {
						output[i * cols + j] = 0;
					}
					break;
				default:
					output[i * cols + j] = 0;
					break;
				}

			}
			else {
				output[i * cols + j] = 0;
			}
		}
	}
}
int* houghTransformLine(unsigned char* binaryimage,int rows,int cols) {
	//uzaklık max köşegen kadar değer alabilir
	int d = (int)sqrt(rows * rows + cols * cols);
	int* houghSpace = new int[d * 180];

	for (int i = 0; i < d * 180; i++) {
		houghSpace[i] = 0;
	}
	//Açıyı üstten aldım
	int distance = 0;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (binaryimage[i * cols + j] == 255) {
				for (int k = 0; k < 180; k++) {
					distance = (i * sin(k * PI / 180) + j * cos(k * PI / 180));
					//cout << "distance: " << distance << endl;
					if (distance >= 0)
						houghSpace[distance * 180 + k] += 1;
				}
			}
		}
	}
	return houghSpace;
}
bool euclideanDistance(int theta, int distance, int* angles, int* distances, int numOfMax) {
	int bias = 30;
	for (int k = 0; k < numOfMax; k++) {

		if (distances[k] != -1 && angles[k] != -1) {
			int epsilon = sqrt(pow((distances[k] - distance), 2) + pow((angles[k] - theta), 2));
			if (epsilon <= bias) {

				return true;
			}
		}
	}
	return false;
}
int* findMaxPoints(int* houghSpace, int numOfMax, int houghHeight, int* angles, int* distances) {
	int* maxPoints = new int[numOfMax];
	for (int i = 0; i < numOfMax; i++) {
		maxPoints[i] = -1;
	}

	int vote = 65;
	int deg;
	int indis = 0;
	bool isSame;
	for (int row = 0; row < houghHeight; row++) {
		for (int col = 0; col < 180; col++) {
			deg = houghSpace[row * 180 + col];
			if (deg > maxPoints[0] && deg > vote) { // values[0] da en küçük max değeri tutuyorum yeni max değerler buldukça kaydırıyorum.
				indis = -1;

				isSame = euclideanDistance(col, row, angles, distances, numOfMax);
				for (int k = 0; k < numOfMax; k++) {
					if (houghSpace[row * 180 + col] > maxPoints[k] && !isSame) { //
						indis = k;
					}
				}
				if (indis != -1) {
					for (int k = 0; k < indis; k++) {
						maxPoints[k] = maxPoints[k + 1];
						angles[k] = angles[k + 1];
						distances[k] = distances[k + 1];
					}
					maxPoints[indis] = houghSpace[row * 180 + col];
					angles[indis] = col;
					distances[indis] = row;
				}
			}
		}
	}


	for (int i = 0; i < numOfMax; i++) {
		cout << "max Point value " << i + 1 << ": " << maxPoints[i] << endl;
		cout << "max Point angle" << i + 1 << ": " << angles[i] << endl;
		cout << "max Point distance " << i + 1 << ": " << distances[i] << endl;
	}
	return maxPoints;
}

void detectLine(unsigned char* image,int* distance, int* angle, int numofline,int cols,int rows) {
	int theta = angle[1];
	int d = distance[1];
	int row = 0;
	int tempdistance = 0;

	for (int r = 0; r < rows; r++) {
		for (int c = 0; c < cols; c++) {
			for (int ang = 0; ang < 180; ang++) {
				tempdistance = (r * sin(ang * PI / 180) + c * cos(ang * PI / 180));
				for (int k = 0; k < numofline; k++) {
					theta = angle[k];
					d = distance[k];
					if (ang == theta && tempdistance == d ) { // && binaryimage[r * cols + c] == 255
						//binaryimage[r * cols + c] = 100;
						image[(r + 1) * (cols + 2) * 3 + (c + 1) * 3 + 2] = 255;
						image[(r + 1) * (cols + 2) * 3 + (c + 1) * 3 + 1] = 0;
						image[(r + 1) * (cols + 2) * 3 + (c + 1) * 3] = 0;
					}
				}
			}
		}
	}
}
int* circleHoughTransform(unsigned char* binaryImage, int radius,int rows,int cols) {
	int* houghSpace = new int[cols * rows];
	for (int i = 0; i < rows * cols; i++) {
		houghSpace[i] = 0;
	} // houghSpace[b][a]
	int a; //a=x0-r*cos(theta)
	int b; //b=y0+r*sin(theta)

	for (int row = 0; row < rows; row++) {
		for (int col = 0; col < cols; col++) {
			if (binaryImage[row * cols + col] == 255) {
				for (int angle = 0; angle < 360;angle++) {
					a = col - radius * cos(angle * PI / 180);
					b = row + radius * sin(angle * PI / 180);
					if (a >= 0 && a < cols && b >= 0 && b < rows)
						houghSpace[b * cols + a]++;
				}
			}
		}
	}
	return houghSpace;
}
int euclidianDistanceforCircle(CircleCenter* centers,int numOfCircle, int x, int y) {
	int bias =10;
	for (int k = 0; k < numOfCircle; k++) {
		if (centers[k].x != 0 && centers[k].y != 0) {
			int epsilon = sqrt(pow((centers[k].x - x), 2) + pow((centers[k].y - y), 2));
			if (epsilon <= bias)
				return true;
		}
	}
	return false;
}
CircleCenter* findMaxValuesforCircle(int* houghSpace, int numOfCircle,int rows,int cols) {
	CircleCenter* centers = new CircleCenter[numOfCircle];
	int* values = new int[numOfCircle];
	for (int i = 0; i < numOfCircle; i++) {
		centers[i].x = 0;
		centers[i].y = 0;
		values[i] = 0;
	}

	int indis = 0;
	int deg = 0;
	int isSame = 0;
	int bias = 30;
	int vote = 80;

	for (int row = 0; row < rows; row++) {
		for (int col = 0; col < cols; col++) {
			deg = houghSpace[row * cols + col];
			if (deg > values[0] && deg > vote) { // values[0] da en küçük max değeri tutuyorum yeni max değerler buldukça kaydırıyorum.
				indis = -1; 
				isSame = euclidianDistanceforCircle(centers, numOfCircle, col, row);
				for (int k = 0; k < numOfCircle; k++) {
					if (houghSpace[row * cols + col] > values[k] && !isSame) { //
						indis = k;
					}
				}
				if (indis != -1) {
					for (int k = 0; k < indis; k++) {
						values[k] = values[k + 1];
						centers[k].x = centers[k + 1].x;
						centers[k].y = centers[k + 1].y;
					}
					values[indis] = houghSpace[row * cols + col];
					centers[indis].x = col;
					centers[indis].y = row;
				}
			}
		}
	}
	for (int i = 0; i < numOfCircle; i++) {
		cout << "values: " << values[i] << endl;
		cout << "x: " << centers[i].x << endl;
		cout << "y: " << centers[i].y << endl;
		cout << endl;
	}
	return centers;
}
void drawSelectedCircle(unsigned char* image, CircleCenter* centers, int numOfCircle, int radius,int rows,int cols) {
	int x;  //x = a+r*cos(theta)
	int y;  //y = b-r*sin(theta)


	for (int k = 0; k < numOfCircle; k++) {
		for (int angle = 0; angle < 360; angle++) {
			x = centers[k].x + radius * cos(angle * PI / 180);
			y = centers[k].y - radius * sin(angle * PI / 180);
			if (x >= 0 && x < cols && y >= 0 && y < rows)
				//binaryImage[y * cols + x] = 100;
				image[(y + 1) * (cols + 2) * 3 + (x + 1) * 3 + 2] = 255;
				image[(y + 1) * (cols + 2) * 3 + (x + 1) * 3 + 1] = 0;
				image[(y + 1) * (cols + 2) * 3 + (x + 1) * 3] = 0;
		}
	}
}
//void say(int* binaryimage, int d) {
//	for (int i = 0; i < d; i++) {
//		for (int j = 0; j < 180;j++) {
//			if (binaryimage[i * 180 + j] > 100) {
//				cout << "Hough degeri: " << binaryimage[i * 180 + j] << endl;
//				cout << "Angle degeri: " << j << endl;
//				cout << "Distance degeri: " << i << endl;
//			}
//		}
//	}
//}
//void say2(unsigned char* binaryimage) {
//	int k = 0;
//	int iki = 0;
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < cols; j++) {
//			if (binaryimage[i * cols + j] != 0 && binaryimage[i * cols + j] != 255 && binaryimage[i * cols + j] != 100) {
//				cout << "görüntü istenilen degil\t i: " << i << "\tj: " << j << "\t" << (int)binaryimage[i * cols + j] << endl;
//				k++;
//			}
//			/*		if (i == 0 && (j == 562 || j == 786 || j == 1089)) {
//						cout << "Degerler\t i: " << i << "\tj: " << j << "\t" << (int)binaryimage[i * cols + j] << endl;
//					}*/
//
//		}
//	}
//	cout << "toplam  olmayan pixel sayısı: " << k << endl;
//	//cout << "toplam binary olan pixel sayısı: " << iki << endl;
//}


//int main()
//{
//	//string image_path = "/Users/elifs/Pictures/Screenshots/line.png";
//	string image_path = "/Users/elifs/Pictures/Screenshots/line3.png";
//	//string image_path = "/Users/elifs/Pictures/Screenshots/resim3.jpg";
//	image = imread(image_path, IMREAD_COLOR);
//	rows = image.rows;
//	cols = image.cols;
//	unsigned char* grayimage = new unsigned char[image.rows * image.cols];
//	RGBtoGray(image.data, grayimage, rows, cols);
//	Mat gray_image(rows, cols, CV_8UC1, grayimage);
//
//	imshow("IMAGE", image);
//	waitKey(0);
//	imshow("Gray Image", gray_image);
//	waitKey(0);
//
//	//int* gradientDirection;
//	int* gradientImage = GradientImage(grayimage,rows,cols);
//	//say(gradientImage,d);
//	unsigned char* normalized_image = NormalizeValue(gradientImage, cols, rows);
//	Mat gradientimage(rows, cols, CV_8UC1, normalized_image);
//	imshow("gradient Image", gradientimage);
//	waitKey(0);
//
//	int* nonmaxSup = new int[rows * cols];
//
//
//	applyNonMaximumSuppression(gradientImage, nonmaxSup,rows,cols);
//
//
//	int* hist = Histogram(nonmaxSup);
//	unsigned char* hysThreshold = new unsigned char[rows * cols];
//	hysteresisThresholding(nonmaxSup, hysThreshold,rows,cols);
//
//	unsigned char* binaryImage = new unsigned char[rows * cols];
//	edgesWithHysteresis(hysThreshold, binaryImage,rows,cols);
//
//	Mat binary_image(rows, cols, CV_8UC1, binaryImage);
//	imshow("BINARY IMAGE", binary_image);
//	waitKey(0);
//	int d = (int)sqrt(rows * rows + cols * cols);
//	int* houghspace = houghTransformLine(binaryImage,rows,cols);
//	//say(houghspace, d);
//	unsigned char* normalizedHough = NormalizeValue(houghspace, 180, d);
//	Mat hough(d, 180, CV_8UC1, normalizedHough);
//	imshow("Hough Sapce", hough);
//	waitKey(0);
//
//	int numOfMax = 10;
//	int* maxPointDistance = new int[numOfMax];
//	int* maxPointAngle = new int[numOfMax];
//	for (int i = 0; i < numOfMax; i++) {
//		maxPointDistance[i] = -1;
//		maxPointAngle[i] = -1;
//	}
//	int* maxPoints = findMaxPoints(houghspace, numOfMax, d, maxPointAngle, maxPointDistance);
//	detectLine(image.data, maxPointDistance, maxPointAngle, numOfMax,cols,rows);
//	
//	imshow("IMAGE WITH LINE", image);
//	waitKey(0);
//}


//unsigned char test[6][6] = { {0,6,4,6,4,2},{2,4,6,5,6,1},{3,4,5,7,4,2},{3,7,6,3,1,2},{2,7,5,5,4,4},{0,7,4,7,0,4} };
//unsigned char ker[3][3] = { {-1,0,1},{-1,0,1} ,{-1,0,1} };
//rows = 6;
//cols = 6;
//Mat testim(6, 6, CV_8UC1, (unsigned char*)test);
//int* result = Convolution(testim);
//int* outputImage = new int[6 * 6];
//for (int i = 0; i < 6 * 6; i++) {
//	outputImage[i] = 0;
//}
//applyNonMaximumSuppression(result, outputImage);
//Circle detection için main fonksiyonu
// 
// 
// 
 
int main() {
	//string image_path = "/Users/elifs/Pictures/Screenshots/circle.jpg";
	string image_path = "/Users/elifs/Pictures/Screenshots/circle2.jpg";
	image = imread(image_path, IMREAD_COLOR);
	rows = image.rows;
	cols = image.cols;
	unsigned char* grayimage = new unsigned char[image.rows * image.cols];
	RGBtoGray(image.data, grayimage, rows, cols);
	Mat gray_image(rows, cols, CV_8UC1, grayimage);

	imshow("IMAGE", image);
	waitKey(0);
	imshow("Gray Image", gray_image);
	waitKey(0);


	int* gradientImage = GradientImage(grayimage,rows,cols);
	//say(gradientImage,d);
	unsigned char* normalized_image = NormalizeValue(gradientImage, cols, rows);
	Mat gradientimage(rows, cols, CV_8UC1, normalized_image);
	imshow("gradient Image", gradientimage);
	waitKey(0);

	int* nonmaxSup = new int[rows * cols];


	applyNonMaximumSuppression(gradientImage, nonmaxSup,rows,cols);


	int* hist = Histogram(nonmaxSup);
	unsigned char* hysThreshold = new unsigned char[rows * cols];
	hysteresisThresholding(nonmaxSup, hysThreshold,rows,cols);

	unsigned char* binaryImage = new unsigned char[rows * cols];
	edgesWithHysteresis(hysThreshold, binaryImage,rows,cols);

	Mat binary_image(rows, cols, CV_8UC1, binaryImage);
	imshow("BINARY IMAGE", binary_image);
	waitKey(0);
	//circle1 78
	int radius = 87;
	int* houghSpace = circleHoughTransform(binaryImage, radius,rows,cols);
	unsigned char* normalizedHough = NormalizeValue(houghSpace, cols, rows);
	Mat hough(rows, cols, CV_8UC1, normalizedHough);
	imshow("hough space", hough);
	waitKey(0);
	int numOfMax = 5;
	CircleCenter* centers = findMaxValuesforCircle(houghSpace, numOfMax,rows,cols);
	drawSelectedCircle(image.data, centers, numOfMax, radius,rows,cols);
	imshow("IMAGE WITH CIRCLE", image);
	waitKey(0);

	/*Canny(image, gray_image, 120, 90);
	imshow("IMAGE", gray_image);
	waitKey(0);*/
	return 0;
}

