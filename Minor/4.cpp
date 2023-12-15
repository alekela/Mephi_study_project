#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


void Convolution1d(double* res, double* data, int len_data, vector<double>& Kernel, int len_Kernel) {
	for (int i = 0; i< len_data + len_Kernel - 1; i++) {
		double sum = 0;
		for (int m = i - len_Kernel + 1; i - m >= 0; m++) {
			sum += data[m] * Kernel[i - m];
		}
		res[i] = sum;
	}
}


void Convolution2D(unsigned char* pSrc, unsigned char* pRes, int iWidth, int iHeight, double* pKernel, int iKernelSize) {
	int iHalf = iKernelSize / 2;
	memcpy(pRes, pSrc, sizeof(*pRes) * iWidth * iHeight);
	// копирование изображения (для полос по краям изображения)
	for (int y = iHalf; y < iHeight - iHalf; ++y) {
		for (int x = iHalf; x < iWidth - iHalf; ++x) {
			double* pk = pKernel;
			const unsigned char* ps = &pSrc[(y - iHalf) * iWidth + x - iHalf];
			int iSum = 0;
			for (int v = 0; v < iKernelSize; ++v) {
				for (int u = 0; u < iKernelSize; ++u)
					iSum += ps[u] * pk[u];

				pk += iKernelSize;  // Переход к следующей строкам
				ps += iWidth;
			}
			iSum = iSum > 0 ? iSum % 256 : 0;
			pRes[y * iWidth + x] = (unsigned char)iSum;
		}
	}
}


void ImageProcessingGray(unsigned char* pOut
	, unsigned char* pIn
	, size_t nWidth
	, size_t nHeight)
{
	const int HPFKernel2d[9] = {1, 1, 1, 
								1, -8, 1, 
								1, 1, 1 };
	int HPFlen = 3;							
	double LPFKernel2d[25] = {1., 4., 7., 4., 1., 
							4., 16., 26., 16., 4., 
							7., 26., 41., 26., 7., 
							4., 16., 26., 16., 4., 
							1., 4., 7., 4., 1.};
	int LPFlen = 5;
	double s = 0.;
	for (int i = 0; i < 25; i++) {
		LPFKernel2d[i] /= 273.;
		s += LPFKernel2d[i];
	}
	std::cout << s;
	Convolution2D(pIn, pOut, nWidth, nHeight, LPFKernel2d, LPFlen);
	return;
}


int main() {
	double test[4] = {2, 1, 3, -1};
	int kernel_len;
	cin >> kernel_len;
	vector <double> kernel(kernel_len);
	for (int i = 0; i < kernel_len; i++) {
		cin >> kernel[i];
	}
	cout << endl;
	double result[6];
	Convolution1d(result, test, 4, kernel, kernel_len);
	for (int i = 0; i < 6; i++) {
		cout << result[i] << endl;
	}



	/*
	double LPFKernel2d[3][3] = {{0.125, 0.125, 0.125}, {0.125, 0, 0.125}, {0.125, 0.125, 0.125}};
	double HPFKernel2d[3][3] = {{-1, -1, -1}, {1, 0, -1}, {1, 1, 1}};
	double LPFKernel1d[3] = {0.5, 0, 0.5};
	double HPFKernel1d[3] = {0.5, 0, -0.5};

	*/
}