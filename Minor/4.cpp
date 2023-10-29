#include <iostream>
#include <vector>

using namespace std;


void Convolution2D(unsigned char* pSrc, unsigned char* pRes, int iWidth, int iHeight, const int* pKernel, int iKernelSize){
    int iHalf = iKernelSize / 2;
    memcpy( pRes, pSrc, sizeof(*pRes) * iWidth * iHeight ); 
	// копирование изображения (для полос по краям изображения)
    for ( int y = iHalf; y < iHeight - iHalf; ++y ) {
		for( int x = iHalf; x < iWidth - iHalf; ++x ){
            const int* pk = pKernel;
            const unsigned char* ps = &pSrc[(y - iHalf) * iWidth + x - iHalf];
            int iSum = 0;
            for ( int v = 0; v < iKernelSize; ++v ){
                for ( int u = 0; u < iKernelSize; ++u )
                    iSum += ps[u] * pk[u];

                pk += iKernelSize;  // Переход к следующей строкам
                ps += iWidth;
            }
			iSum = iSum > 0 ? iSum % 256 : 0;
            pRes[y * iWidth + x] = (unsigned char)iSum;
        }
    }
}


void Convolution1d(double* res, double* data, int len_data, vector<double>& Kernel, int len_Kernel) {
	for (int i = 0; i< len_data + len_Kernel - 1; i++) {
		double sum = 0;
		for (int m = 0; m < len_data; m++) {
			if (i - m >= 0 && i - m < len_Kernel) {
				sum += data[m] * Kernel[i - m];
			}
		}
		res[i] = sum;
	}
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