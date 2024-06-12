// gorunutIsleme.cpp : Bu dosya 'main' işlevi içeriyor. Program yürütme orada başlayıp biter.
//
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <vector>
#include "Object.h"
#include <stdlib.h>
#include "sqlite3.h"
#include "myHuMoments.h"
#include <math.h>

using namespace std;
unsigned char* input;
struct Nesne {
	int rowStart;
	int rowFinish;
	int colStart;
	int colFinish;
	unsigned char* object;

	Nesne(int rows, int rowf, int cols, int colf, unsigned char* obj) :rowStart(rows), rowFinish(rowf), colStart(cols), colFinish(colf), object(obj) {}
};

void RGBtoGray(unsigned char* input, unsigned char* grayimage, int rows, int cols) {

	for (int r = 0; r < rows; r++) {
		for (int c = 0; c < cols; c++) {
			grayimage[r * cols + c] = input[r * cols * 3 + c * 3 + 2] * 0.299 + input[r * cols * 3 + c * 3 + 1] * 0.587 + input[r * cols * 3 + c * 3] * 0.114;
		}
	}
}
int* Histogram(unsigned char* inp, int h, int w) {
	int* hist = new int[256];
	for (int i = 0; i < 256; i++) {
		hist[i] = 0;
	}
	for (int i = 0; i < h * w; i++) {
		hist[inp[i]]++;
	}
	//delete[] hist;
	return hist;
}

unsigned char* kMeans(unsigned char* inp, int r, int c, int* hist, float k1, float k2) {
	unsigned char* out = new unsigned char[r * c];
	float K11 = 0;
	float K22 = 0;
	int counter1 = 0;  //k1'e yakın kaç pixel olduğunu tutar
	int counter2 = 0;  //k2'ye yakın kaç pixel olduğunu tutar



	for (int i = 0; i < 256; i++) {
		if (abs(k1 - i) < abs(k2 - i)) {
			K11 += i * hist[i];
			counter1 += hist[i];
		}
		else {
			K22 += i * hist[i];
			counter2 += hist[i];
		}
	}
	K11 = K11 / counter1;
	K22 = K22 / counter2;

	if (K11 == k1 && K22 == k2) {
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				if (abs(inp[i * c + j] - K11) < abs(inp[i * c + j] - K22)) {
					out[i * c + j] = 0;  // Değeri k2 ye yakın olanlara siyah atandı bunlar zemin
				}
				else {
					out[i * c + j] = 255;  // Değeri k1 ye yakın olanlara beyaz atandı bunlar nesne
					
				}
			}
		}
		return out;

	}
	else {
		return kMeans(inp, r, c, hist, K11, K22);
	}
}

unsigned char* Dilation(unsigned char* inp, int r, int c) {
	unsigned char kernel[3][3] = { {0,255,0},{255,255,255},{0,255,0} };
	unsigned char* structurElement = (unsigned char*)kernel;
	//unsigned char* output = new unsigned char[r*c];
	unsigned char* output = new unsigned char[r * c];

	//unsigned char* ninput = new unsigned char[(r + 2) * (c + 2)];

	//for (int i = 0; i < (r + 2); i++) {
	//	for (int j = 0; j < (c + 2); j++) {
	//		if (i == 0 || i == r + 1 || j == 0 || j == c + 1) {
	//			ninput[i * (c + 2) + j] = 255;
	//		}
	//		else {
	//			ninput[i * (c + 2) + j] = inp[(i - 1) * c + j - 1];
	//		}
	//	}
	//}
	//
	//for (int i = 0; i < r; i++) {
	//	for (int j = 0; j < c; j++) {
	//		output[i * c + j] = 255;
	//		for (int m = 0; m < 3; m++) {
	//			for (int n = 0; n < 3; n++) {
	//				
	//				int result = ninput[(i + m) * c + n];
	//				if (result == 0 && ninput[(i + m) * c + n] == structurElement[m * 3 + n]) {
	//					output[i * c + j] = 0;
	//				}
	//			}
	//		}
	//	}
	//}

	//return output;
	////return ninput;

	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			output[i * c + j] = 0; // Çıktı görüntüsünü sıfırla
			for (int m = -1; m <= 1; m++) {
				for (int n = -1; n <= 1; n++) {
					if (i + m >= 0 && i + m < r && j + n >= 0 && j + n < c && structurElement[(m + 1) * 3 + n + 1] == 0) {
						// Yapılandırma elemanı içindeki pikselleri kontrol et ve dilatasyon işlemi uygula
						output[i * c + j] = std::max(output[i * c + j], inp[(i + m) * c + j + n]);
					}
				}
			}
		}
	}
	//erosion olarak çalıştı bende şaşkınım!!
	return output;
}

vector<int> Labeling(unsigned char* inp, int r, int c, unsigned char* out) {

	unsigned char* ninput = new unsigned char[(r + 2) * (c + 2)];
	int neighbours[4];
	//int tag = 2;
	int tag = 200;
	int count = 0;
	vector<int> tags;  //Çerçeve belirlemek için etiketlerin ne olduğunu tutar.

	//işlem kolaylığı için padding verdim
	for (int i = 0; i < (r + 2); i++) {
		for (int j = 0; j < (c + 2); j++) {
			if (i == 0 || i == r + 1 || j == 0 || j == c + 1) {
				ninput[i * (c + 2) + j] = 0;
			}
			else {
				ninput[i * (c + 2) + j] = inp[(i - 1) * c + j - 1];
			}
		}
	}
	//Birinci aşama
	for (int i = 1; i < r + 1; i++) {
		for (int j = 1; j < c + 1; j++) {
			int nzero = 0;
			if (ninput[i * (c + 2) + j] == 255) {
				neighbours[0] = ninput[(i - 1) * (c + 2) + j - 1];
				neighbours[1] = ninput[(i - 1) * (c + 2) + j];
				neighbours[2] = ninput[(i - 1) * (c + 2) + j + 1];
				neighbours[3] = ninput[i * (c + 2) + j - 1];

				for (int k = 0; k < 4; k++) {
					if (neighbours[k] != 0) {
						nzero++;
					}
				}
				if (nzero == 0) {
					ninput[i * (c + 2) + j] = tag;
					tags.push_back(tag);
					//cout << count << ". Etiket: " << tag << endl;
					tag -= 2;
					//tag++;
					count++;
				}
				else if (nzero == 1) {
					for (int k = 0; k < 4; k++) {
						if (neighbours[k] != 0) {
							ninput[i * (c + 2) + j] = neighbours[k];
							break;
						}
					}
				}
				else {
					//int min = 256;
					int max = 0;
					int change = -1;
					for (int k = 0; k < 4; k++) {
						if (neighbours[k] != 0 && neighbours[k] > max) {
							max = neighbours[k];
						}
					}
					ninput[i * (c + 2) + j] = max;
					for (int k = 0; k < 4; k++) {
						if (neighbours[k] != 0 && neighbours[k] != max) {
							change = neighbours[k];
							break;
						}
					}
					if (change != -1) {
						count--;
						cout << "Silinen Etiket :" << change << endl;

						//vektöre eklenen change değerli etiketi silme
						for (int k = 0; k < tags.size(); k++) {
							if (tags[k] == change) {
								tags.erase(tags.begin() + k);
								break;
							}

						}

						//matris üzerinde change'e eşit tüm değerler max değer yapılır.
						for (int m = 1; m <= i; m++) {
							for (int n = 1; n <= c; n++) {
								if (ninput[m * (c + 2) + n] == change) {
									ninput[m * (c + 2) + n] = max;
								}
							}
						}
					}

				}
			}
		}
	}
	//görüntü üzerindeki paddingi kaldırma 
	for (int i = 1; i < (r + 1); i++) {
		for (int j = 1; j < (c + 1); j++) {
			out[(i - 1) * c + j - 1] = ninput[i * (c + 2) + j];
		}
	}

	return tags;
}

double* Moments1(unsigned char* object, int row, int col) {
	double* moments = new double[7];
	double M00 = 0.0;
	double M10 = 0.0;
	double M01 = 0.0;
	for (int x = 0; x < col; x++) {
		for (int y = 0; y < row; y++) {
			M00 += pow(x, 0) * pow(y, 0) * object[y * col + x];
			M10 += pow(x, 1) * pow(y, 0) * object[y * col + x];
			M01 += pow(x, 0) * pow(y, 1) * object[y * col + x];
		}
	}
	double centerX = M10 / M00;
	double centerY = M01 / M00;

	//Merkezi momentler U ile temsil edildi.
	double U00 = 0;
	double U20 = 0;
	double U02 = 0;
	double U11 = 0;
	double U30 = 0;
	double U12 = 0;
	double U21 = 0;
	double U03 = 0;

	for (int x = 0; x < col; x++) {
		for (int y = 0; y < row; y++) {
			U00 += pow(x - centerX, 0) * pow(y - centerY, 0) * object[y * col + x];
			U20 += pow(x - centerX, 2) * pow(y - centerY, 0) * object[y * col + x];
			U02 += pow(x - centerX, 0) * pow(y - centerY, 2) * object[y * col + x];
			U11 += pow(x - centerX, 1) * pow(y - centerY, 1) * object[y * col + x];
			U30 += pow(x - centerX, 3) * pow(y - centerY, 0) * object[y * col + x];
			U12 += pow(x - centerX, 1) * pow(y - centerY, 2) * object[y * col + x];
			U21 += pow(x - centerX, 2) * pow(y - centerY, 1) * object[y * col + x];
			U03 += pow(x - centerX, 0) * pow(y - centerY, 3) * object[y * col + x];
		}
	}
	//Normalize merkezi moment N ile temsil edildi.
	double N20 = U20 / pow(U00, 2);
	double N02 = U02 / pow(U00, 2);
	double N11 = U11 / pow(U00, 2);
	double N30 = U30 / pow(U00, (3 / 2) + 1);
	double N12 = U12 / pow(U00, (3 / 2) + 1);
	double N21 = U21 / pow(U00, (3 / 2) + 1);
	double N03 = U03 / pow(U00, (3 / 2) + 1);


	//yedi değişmmez moment kümesi
	moments[0] = N20 + N02;
	moments[1] = pow((N20 - N02), 2) + 4 * pow(N11, 2);
	moments[2] = pow((N30 - 3 * N12), 2) + pow((3 * N21 - N03), 2);
	moments[3] = pow((N30 + N12), 2) + pow((N21 + N03), 2);
	moments[4] = (N30 - 3 * N12) * (N30 + N12) * (pow((N30 + N12), 2) - 3 * pow((N21 + N03), 2)) + (3 * N21 - N03) * (N21 + N03) * (3 * pow((N30 + N12), 2) - pow(N21 + N03, 2));
	moments[5] = (N20 - N02) * (pow((N30 + N12), 2) - pow(N21 + N03, 2)) + 4 * N11 * (N30 + N12) * (N21 + N03);
	moments[6] = (3 * N21 - N03) * (N30 + N12) * (pow((N30 + N12), 2) - 3 * pow((N21 + N03), 2)) + (3 * N12 - N30) * (N21 + N03) * (3 * pow((N30 + N12), 2) - pow((N21 + N03), 2));

	/*cout << "m00: " << M00 << endl;
cout << "m10: " << M10 << endl;
cout << "m01: " << M01 << endl;
cout << "r: " << centerX << endl;
cout << "c: " << centerY << endl;
cout << "F1: " << moments[0] << endl;
cout << "F2: " << moments[1] << endl;
cout << "F3: " << moments[2] << endl;
cout << "F4: " << moments[3] << endl;
cout << "F5: " << moments[4] << endl;
cout << "F6: " << moments[5] << endl;
cout << "F7: " << moments[6] << endl;*/
	
	return moments;
}


Nesne Framing(unsigned char* input,unsigned char* labeledinp,unsigned char* binaryinp, int height, int width, int tag) {
	int rows = height;
	int rowf = 0;
	int cols = width;
	int colf = 0;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (labeledinp[i * width + j] == tag) {
				if (i < rows) {
					rows = i;
				}
				if (i > rowf) {
					rowf = i;
				}
				if (j < cols) {
					cols = j;
				}
				if (j > colf) {
					colf = j;
				}
			}
		}
	}
	unsigned char* object = new unsigned char[(rowf - rows + 1) * (colf - cols + 1)];

	for (int i = rows; i <= rowf; i++) {
		for (int j = cols; j <= colf; j++) { 
			object[(i - rows) * (colf - cols + 1) + j - cols] = binaryinp[i * width + j];
			if (i == rows || i == rowf || j == cols || j == colf) {
				input[i * width * 3 + j * 3 + 2] = 198;			  
				input[i * width * 3 + j * 3 + 1] = 226;			  
				input[i * width * 3 + j * 3] = 255;			  
			}

		}
	}
	struct Nesne obj = Nesne(rows, rowf, cols, colf, object);
	return	obj;				
}

double calculateDistance(double* object, double* ref) {
	double sum = 0;
	double temp;

	for (int i = 0; i < 7; i++) {
		temp = object[i] - ref[i];
		sum += pow(temp, 2);
	}
	return sqrt(sum);
}

void Classification(int red, int green, int blue,Nesne nesne,int width,unsigned char* copyImage) {

	for (int i = nesne.rowStart; i <= nesne.rowFinish; i++) {
		for (int j = nesne.colStart; j <= nesne.colFinish; j++) {
			if(nesne.object[(i - nesne.rowStart) * (nesne.colFinish - nesne.colStart + 1) + j - nesne.colStart] == 255){
			copyImage[i * width * 3 + j * 3 + 2] = red;
			copyImage[i * width * 3 + j * 3 + 1] = green;
			copyImage[i * width * 3 + j * 3] = blue;
			}

		}
	}
}
unsigned char* CopyImage(unsigned char* image, int width, int height) {
	unsigned char* copyIm = new unsigned char[width * height * 3];

	for (int r = 0; r < height*width*3; r++) {
		copyIm[r] = image[r];
	}
	return copyIm;
}
double detectMin(double* distances, int numOfRef) {
	double min = distances[0];

	for (int i = 0; i < numOfRef; i++) {
		if (distances[i] < min)
			min = distances[i];
	}
	return min;
}
sqlite3* db;
sqlite3_stmt* stmt;
int result;

void connection();
void insertData();
void readData(string objectName, double* moment);
using namespace cv;
using namespace std;


int main()
{
	double* kareMoment = new double[7];
	double* ucgenMoment = new double[7];
	double* artiMoment = new double[7];
	int kare = 0;
	int ucgen = 0;
	int arti = 0;
	
	connection();
	//insertData();
	readData("kare", kareMoment);
	readData("ucgen", ucgenMoment);
	readData("arti", artiMoment);
	cout << "kare" << endl;
	for (int i = 0; i < 7; i++) {
		cout << i << ": " << kareMoment[i] << endl;
	}
	cout << endl;
	cout << "ucgen" << endl;
	for (int i = 0; i < 7; i++) {
		cout << i << ": " << ucgenMoment[i] << endl;
	}
	cout << endl;
	cout << "arti" << endl;
	for (int i = 0; i < 7; i++) {
		cout << i << ": " << artiMoment[i] << endl;
	}

	//string image_path = "/Users/elifs/Pictures/Screenshots/resim64.jpg";
	string image_path = "/Users/elifs/Pictures/Screenshots/testImage2.jpg";
	Mat image = imread(image_path, IMREAD_COLOR);
	input = image.data;
	int rows = image.rows;
	int cols = image.cols;
	//Mat imagecopy = image.clone();
	unsigned char* copyImage = CopyImage(input, cols, rows);
	Mat copy_image(rows,cols, CV_8UC3,copyImage);

	imshow("IMAGE", image);
	waitKey(0);

	

	unsigned char* grayimage = new unsigned char[image.rows * image.cols];
	RGBtoGray(input, grayimage, rows, cols);
	Mat gray_image(rows, cols, CV_8UC1, grayimage);
	imshow("Gray Image", gray_image);
	waitKey(0);

	int* hist = Histogram(grayimage, rows, cols);

	unsigned char* binaryimage = kMeans(grayimage, rows, cols, hist, 0, 200);
	Mat binary_image(rows, cols, CV_8UC1, binaryimage);
	imshow("Binary Image", binary_image);
	waitKey(0);

	unsigned char* label = new unsigned char[image.rows * image.cols];
	vector<int> tags = Labeling(binaryimage, rows, cols, label);
	int red = 198;
	int green = 226;
	int blue = 255;
	double* distances = new double[3];
	for (int k = 0; k < tags.size(); k++) {
		Nesne obj = Framing(input,label,binaryimage, rows, cols, tags[k]);  
		double* objectmoment = Moments1(obj.object, obj.rowFinish - obj.rowStart + 1, obj.colFinish - obj.colStart + 1); // bilinmeyen nesnenin moment değeri
		distances[0] = calculateDistance(objectmoment, kareMoment);
		distances[1] = calculateDistance(objectmoment, ucgenMoment);
		distances[2] = calculateDistance(objectmoment, artiMoment);
		cout << k+1<<". resim:" << endl;
		cout << "Ucgen distance: " << distances[1] << endl;
		cout << "Kare distance: " << distances[0] << endl;
		cout << "Arti distance: " << distances[2] << endl;
		 
		double min = detectMin(distances, 3);
		if (distances[0] == min) {
			cout << "Cisim karedir" << endl;
			kare++;
			red = 255;
			green = 0;
			blue = 0;
		}
		else if (distances[1] == min) {
			cout << "Cisim ucgendir" << endl;
			ucgen++;
			red = 0;
			green = 255;
			blue = 0;
		}
		else {
			cout << "artıdır " << endl;
			arti++;
			red = 0;
			green = 0;
			blue = 255;
		}
		Classification(red, green, blue, obj, cols,copyImage);
		cout << endl;
		cout << endl;
		cout << endl;
		cout << endl;
		cout << endl;
		
		
	}
	int temp = 0;
	if (arti < 3) temp += 3 - arti;
	if (kare < 4) temp += 4 - kare;
	if (ucgen < 3) temp += 3 - ucgen;
	cout << "Basari orani: %" << (10 - temp) * 10;
	imshow("IMAGE AFTER FRAMING", image);
	waitKey(0);
	imshow("IMAGE AFTER CLASSIFICATION Yeşil(üçgen) - Kırmızı(kare) - Mavi(artı)", copy_image);
	waitKey(0);
	delete[] grayimage;
	delete[] label;


	return 0;
}
void connection() {
	if (sqlite3_open("moments.db", &db) == SQLITE_OK) {
		result = sqlite3_prepare_v2(db, " CREATE TABLE IF NOT EXISTS moments (objectName varchar(80),F1	double, F2 double,F3	double, F4	double,F5 double,F6	double, F7 double);", -1, &stmt, NULL);
		sqlite3_step(stmt);
		sqlite3_finalize(stmt);

		if (result != SQLITE_OK) {
			cout << "Error: " << sqlite3_errmsg(db) << "\n";
		}
		else {
			cout << "Table created successfully";
		}
	}
}
void insertData() {
	string image_path = "/Users/elifs/Pictures/Screenshots/samplearti.jpg";
	double momentlerdb[7] = { 0.0 };
	//readData("mercimek", momentlerdb);
	Mat image = imread(image_path, IMREAD_COLOR);
	input = image.data;
	unsigned char* grayimage = new unsigned char[image.rows * image.cols];
	int rows = image.rows;
	int cols = image.cols;

	imshow("IMAGE", image);
	waitKey(0);


	RGBtoGray(input, grayimage, rows, cols);
	Mat gray_image(rows, cols, CV_8UC1, grayimage);
	namedWindow("Gray Image", WINDOW_AUTOSIZE);
	imshow("Gray Image", gray_image);
	waitKey(0);

	int* hist = Histogram(grayimage, rows, cols);

	unsigned char* binaryimage = kMeans(grayimage, rows, cols, hist, 0, 200);
	Mat binary_image(rows, cols, CV_8UC1, binaryimage);

	namedWindow("Binary Image", WINDOW_AUTOSIZE);
	imshow("Binary Image", binary_image);
	waitKey(0);


	/*Mat kernel = getStructuringElement(MORPH_RECT, Size(5, 5));

	dilate(binary_image, binary_image, kernel);
	namedWindow("Dilated Image", WINDOW_NORMAL);
	imshow("dilate", binary_image);
	waitKey(0);*/
	unsigned char* label = new unsigned char[image.rows * image.cols];
	vector<int> tags = Labeling(binaryimage, rows, cols, label);
	int r = 198;
	int g = 226;
	int b = 255;

	cout << "nesne sayısı: " << tags.size() << endl;
	for (int k = 0; k < tags.size(); k++) {
		Nesne obj = Framing(input,label, binaryimage, rows, cols, tags[k]);  // bilinmeyen nesnenin moment değeri
		double* objectmoment = Moments1(obj.object, obj.rowFinish - obj.rowStart + 1, obj.colFinish - obj.colStart + 1);

		for (int i = 0; i < 7; i++) {
			momentlerdb[i] += objectmoment[i];
		}
		
	}
	for (int i = 0; i < 7; i++) {
		momentlerdb[i] = momentlerdb[i] / 5;
	}
	string objectname = "arti";
	string query = "INSERT INTO moments(objectName,F1,F2,F3,F4,F5,F6,F7) VALUES(?,?,?,?,?,?,?,?);";
	result = sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, NULL);
	sqlite3_bind_text(stmt, 1, objectname.c_str(), objectname.length(), SQLITE_TRANSIENT);
	sqlite3_bind_double(stmt, 2, momentlerdb[0]);
	sqlite3_bind_double(stmt, 3, momentlerdb[1]);
	sqlite3_bind_double(stmt, 4, momentlerdb[2]);
	sqlite3_bind_double(stmt, 5, momentlerdb[3]);
	sqlite3_bind_double(stmt, 6, momentlerdb[4]);
	sqlite3_bind_double(stmt, 7, momentlerdb[5]);
	sqlite3_bind_double(stmt, 8, momentlerdb[6]);

	sqlite3_step(stmt);
	sqlite3_finalize(stmt);

	if (result != SQLITE_OK) {
		cout << "Error: " << sqlite3_errmsg(db) << "\n";
	}
	else {
		cout << "Table inserted successfully";
	}
	imshow("Training Objects", image);
	waitKey(0);

	delete[] grayimage;
	delete[] label;
}
void readData(string objectName, double* moment) {

	string query = "Select * from moments WHERE objectName = ?";
	result = sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, NULL);
	sqlite3_bind_text(stmt, 1, objectName.c_str(), objectName.length(), SQLITE_TRANSIENT);
	if (result != SQLITE_OK) {
		cout << "Error: " << sqlite3_errmsg(db) << "\n";
	}
	else {
		while ((result = sqlite3_step(stmt)) == SQLITE_ROW) {
			moment[0] = sqlite3_column_double(stmt, 1);
			moment[1] = sqlite3_column_double(stmt, 2);
			moment[2] = sqlite3_column_double(stmt, 3);
			moment[3] = sqlite3_column_double(stmt, 4);
			moment[4] = sqlite3_column_double(stmt, 5);
			moment[5] = sqlite3_column_double(stmt, 6);
			moment[6] = sqlite3_column_double(stmt, 7);

		}
		
	}
}

