#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <complex.h>
#include <vector>


using namespace std;


void print(double* data, int len) {
    for (int i = 0; i < len; i++) {
        cout << data[i] << " ";
    }
    cout << "\n";
}


void quantization(double* data, int len_data, double* quants, int levels) {
    for (int j = 0; j < len_data; j++) {
        double up, down;
        for (int i = 1; i < levels + 1; i++){
            if (data[j] <= quants[i]){
                down = quants[i - 1];
                up = quants[i];
                break;
            }
        }
        data[j] = data[j] - down < up - data[j] ? down : up;
    }
}


double interpolate0(double* data, double* grid, int len_data, double x) {
    if (x < grid[0]) {
        return data[0]; // если точка левее интервала, то возвращается крайнее левое значение функции
    }
    for (int i = 1; i < len_data; i++) {
        if (grid[i] >= x) {
            int k = abs(grid[i] - x) > abs(grid[i - 1] - x) ? i - 1 : i;
            return data[k];
        }
    }
    return data[len_data - 1]; // если точка правее интервала, то возвращается крайнее правое значение функции
}


double interpolate1(double* data, double* grid, int len_data, double x) {
    if (x < grid[0]) {
        return data[0]; // если точка левее интервала, то возвращается крайнее левое значение функции
    }
    for (int i = 1; i < len_data; i++) {
        if (grid[i] >= x) {
            return data[i - 1] + (data[i] - data[i - 1]) / (grid[i] - grid[i - 1]) * (x - grid[i - 1]);
        }
    }
    return data[len_data - 1]; // если точка правее интервала, то возвращается крайнее правое значение функции
}


double my_sin(double x, double period, double phi){
    return sin(x * 2 * 3.1415 / period + phi);
}

int main() {
    double left = 0;
    double right = 10;
    const int N_samples = 5;
    double h = (right - left) / (N_samples - 1);
    double data[N_samples];
    double x_grid[N_samples];
    for (int i = 0; i < N_samples; i++) {
        data[i] = my_sin(left + i * h, 4, 0);
        x_grid[i] = left + i * h;
    }
    const int levels = 5;
    double quants[levels];
    for (int i = 0; i < levels; i++) {
        quants[i] = -1 + 2. / (levels - 1) * i;
    }
    quantization(data, N_samples, quants, levels);

    const int N2_samples = N_samples * 2 - 1;
    double interp_data[N2_samples];
    double h2 = (right - left) / (N_samples - 1) / 2.;
    for (int i = 0; i < N2_samples; i++) {
        interp_data[i] = interpolate0(data, x_grid, N_samples, left + i * h2);
    }
    print(data, N_samples);
    print(interp_data, N2_samples);
}



/*
int main(int argc, char* argv[])
{
	class CBitsPtrGuard
	{
	public:
		CBitsPtrGuard(unsigned char** pB) : m_ppBits(pB) { }
		~CBitsPtrGuard() { if (*m_ppBits) delete* m_ppBits, * m_ppBits = 0; }
	protected:
		unsigned char** m_ppBits;
	};

	// parse input parameters
	char	szInputFileName[256];
	char    szOutputFileName[256];
	if (argc < 2)
		printf("\nformat: pngtest <input_file> [<output_file>]");
	else
	{
		strcpy(szInputFileName, argv[1]);
		if (argc > 2)
			strcpy(szOutputFileName, argv[2]);
		else
		{
			strcpy(szOutputFileName, szInputFileName);
			strcat(szOutputFileName, "_out.png");
		}
	}


	size_t nReqSize = NPngProc::readPngFile(szInputFileName, 0, 0, 0, 0);
	if (nReqSize == NPngProc::PNG_ERROR)
	{
		printf("\nError ocured while pngfile was read");
		return -1;
	}


	unsigned char* pInputBits = new unsigned char[nReqSize];
	if (!pInputBits)
	{
		printf("\nCan't allocate memory for image, required size is %u", nReqSize);
		return -1;
	}
	CBitsPtrGuard InputBitsPtrGuard(&pInputBits);


	unsigned char* pOutputBits = new unsigned char[4 * nReqSize];
	if (!pOutputBits)
	{
		printf("\nCan't allocate memory for image, required size is %u", 4 * nReqSize);
		return -1;
	}


	CBitsPtrGuard OutputBitsPtrGuard(&pOutputBits);

	size_t nWidth, nHeight;
	unsigned int nBPP;

	size_t nRetSize = NPngProc::readPngFileGray(szInputFileName, pInputBits, &nWidth, &nHeight); // last argument may be &nBPP
	nBPP = 8;


	// ASSERT(nRetSize == nReqSize);

	// TODO: image processing 
	ImageProcessingGray(pOutputBits, pInputBits, nWidth, nHeight);

	if (NPngProc::writePngFile(szOutputFileName, pOutputBits, 2 * nWidth, 2 * nHeight, nBPP) == NPngProc::PNG_ERROR)
	{
		printf("\nError ocuured during png file was written");
		return -1;
	}

	return 0;
}


void ResizeBilinear(const unsigned char* pbIn, int lWidthIn,
	int lHeightIn, unsigned char* pbOut, int lWidthOut, int lHeightOut) {
	for (int i = 0; i < lHeightOut; ++i) {
		double yy = (double)i * (double)lHeightIn / lHeightOut;
		int y = (int)yy; // целая часть yy
		double u = yy - (double)y; // дробная часть yy
		for (int j = 0; j < lWidthOut; ++j) {
			double xx = (double)j * (double)lWidthIn / lWidthOut;
			int x = (int)xx; // целая часть xx
			double v = xx - (double)x; // дробная часть xx
			// значения соседних пикселов
			int lP[2][2]{ pbIn[y * lWidthIn + x], pbIn[y * lWidthIn + x + 1],
			pbIn[(y + 1) * lWidthIn + x], pbIn[(y + 1) * lWidthIn + x + 1] };
			// билинейная интерполяция
			pbOut[i * lWidthOut + j] = (unsigned char)((1. - u) * (1. - v) * lP[0][0] +
				u * (1. - v) * lP[1][0] + v * (1. - u) * lP[0][1] +
				u * v * lP[1][1]);
		}
	}
}


void ResizeNN(const unsigned char* pbIn, int lWidthIn,
	int lHeightIn, unsigned char* pbOut, int lWidthOut, int lHeightOut) {
	for (int i = 0; i < lHeightOut; ++i) {
		double yy = (double)i * (double)lHeightIn / lHeightOut;
		int y = (int)yy; // целая часть yy
		double u = yy - (double)y; // дробная часть yy
		for (int j = 0; j < lWidthOut; ++j) {
			double xx = (double)j * (double)lWidthIn / lWidthOut;
			int x = (int)xx; // целая часть xx
			double v = xx - (double)x; // дробная часть xx
			// значения соседних пикселов
			int lP[2][2]{ pbIn[y * lWidthIn + x], pbIn[y * lWidthIn + x + 1],
			pbIn[(y + 1) * lWidthIn + x], pbIn[(y + 1) * lWidthIn + x + 1] };
			// интерполяция ближайшего соседа
			pbOut[i * lWidthOut + j] = (unsigned char)(lP[(int)(xx + 0.5) - x][(int)(yy + 0.5) - y]);
		}
	}
}


void ImageProcessingGray(unsigned char* pOut
	, unsigned char* pIn
	, size_t nWidth
	, size_t nHeight)
{
	ResizeNN(pIn, nWidth, nHeight, pOut, 2 * nWidth, 2 * nHeight);
	return;
}*/