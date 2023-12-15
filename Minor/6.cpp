#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <ctime>
#include <sstream>
#include <fstream>

using namespace std;


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


double my_sin(double x, double period, double phi){
    return sin(x * 2 * 3.1415 / period + phi);
}


double median1d(vector<double> tmp, int len) {
    for (int i = 0; i < len; i++) {
        double t;
        for (int j = 0; j < len - i - 1; j++) {
            if (tmp[j] > tmp[j + 1]) {
                t = tmp[j];
                tmp[j] = tmp[j + 1];
                tmp[j + 1] = t;
            }
        }
    }
    if (len % 2) {
        return tmp[len / 2];
    }
    else {
        return (tmp[len / 2] + tmp[len / 2 - 1]) / 2.;
    }
}


void filter1d(double* res, double* data, int len_data, vector<double> aperture) {
    vector<double> tmp(aperture.size());
    memcpy(res, data, sizeof(*res) * len_data);
    for (int i = 0; i < len_data - aperture.size() + 1; i++) {
        for (int j = 0; j < aperture.size(); j++) {
            tmp[j] = aperture[j] * data[i + j];
        }
        res[i] = median1d(tmp, aperture.size());
    } 
}


void add_impulse_noise(double P, double* data, int len_data, double right, double left) {
    double check;
    double a = right;
    srand(time(NULL));
    for (int i = 0; i < len_data; i++) {
        check = rand() / (double) RAND_MAX;
        if (check < P) {
            if (check < P / 2.) {
                data[i] = right;
            }
            else{
                data[i] = left;
            }
        }
    }
}


int main() {
    // создание отквантованного одномерного сигнала
    double left = 0;
    double right = 10;
    const int N_samples = 500;
    double h = (right - left) / (N_samples - 1);
    double data[N_samples];
    double x_grid[N_samples];
    for (int i = 0; i < N_samples; i++) {
        data[i] = my_sin(left + i * h, 4, 0);
        x_grid[i] = left + i * h;
    }
    const int levels = 100;
    double quants[levels];
    for (int i = 0; i < levels; i++) {
        quants[i] = -1 + 2. / (levels - 1) * i;
    }
    quantization(data, N_samples, quants, levels);
    add_impulse_noise(0.3, data, N_samples, quants[levels - 1], quants[0]);

    const char csv_file_name[64] = "work6.csv";
    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal\n";
    for (size_t i = 0; i < N_samples; ++i){
        csv_file << x_grid[i] << "," << data[i] << "\n";
    }
    csv_file.close();

    int len_ap;
    cin >> len_ap;
    vector<double> aperture;
    double tmp;
    for (int i = 0; i < len_ap; i++) {
        cin >> tmp;
        aperture.push_back(tmp);
    }

    double result[N_samples];
    filter1d(result, data, N_samples, aperture);

    const char csv_file_name2[64] = "work62.csv";
    csv_file.open(csv_file_name2);
    csv_file << "time,signal\n";
    for (size_t i = 0; i < N_samples; ++i){
        csv_file << x_grid[i] << "," << result[i] << "\n";
    }
    csv_file.close();
}  


int median2d(int* tmp, int len) {
	for (int i = 0; i < len; i++) {
		int t;
		for (int j = 0; j < len - i - 1; j++) {
			if (tmp[j] > tmp[j + 1]) {
				t = tmp[j];
				tmp[j] = tmp[j + 1];
				tmp[j + 1] = t;
			}
		}
	}
	if (len % 2) {
		return tmp[len / 2];
	}
	else {
		return (tmp[len / 2] + tmp[len / 2 - 1]) / 2;
	}
}

void ImageProcessingGray(unsigned char* pOut
	, unsigned char* pIn
	, size_t nWidth
	, size_t nHeight)
{
	int aperture_size = 4;
	int aperture[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	int iHalf = aperture_size / 2;
	memcpy(pOut, pIn, sizeof(*pOut) * nWidth * nHeight);
	// копирование изображения (для полос по краям изображения)
	for (int y = iHalf; y < nHeight - iHalf; ++y) {
		for (int x = iHalf; x < nWidth - iHalf; ++x) {
			int pk[16]{};
			const unsigned char* ps = &pIn[(y - iHalf) * nWidth + x - iHalf];
			for (int v = 0; v < aperture_size; ++v) {
				for (int u = 0; u < aperture_size; ++u) {
					pk[u + v * aperture_size] = ps[u] * aperture[u + v * aperture_size];
				}
				ps += nWidth;
			}
			pOut[y * nWidth + x] = (unsigned char) median2d(pk, aperture_size * aperture_size);
		}
	}
	return;
}