#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <random>
#include "windows.h"
#include <math.h>
#include <stdio.h>


using namespace std;


void quantization(double* data, int len_data, double* quants, int levels){

    double up, down;
    for (int j = 0; j < len_data; j++) {
        for (int i = 1; i < levels; i++){
            if (data[j] <= quants[i]){
                down = quants[i - 1];
                up = quants[i];
                break;
            }
        }
        if (data[j] - down < up - data[j]) {
            data[j] = down;
        }
        else {
            data[j] = up;
        }
    }
}


double my_sin(double x, double period, double phi){
    return sin(x * 2 * 3.1415 / period + phi);
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


double gauss(double x, double sigma) {
    return 1. / sqrt(2 * 3.141592) / sigma * exp(- (x * x / 2. / sigma / sigma));
}


void add_Gauss_noise(double disp, double* data, int len_data, double* quants, int levels) {
    double check;
    const int n = 1000;
    double laplas[n];
    double laplas_x[n];
    double h = 10. * disp / (n - 1);
    for (int i = 0; i < n; i++) {
        laplas_x[i] = -5. * disp + i * h;
        if (i == 0) {
            laplas[i] = 0;
        }
        else {
        laplas[i] = laplas[i - 1] + gauss(laplas_x[i], disp) * h;
        }
    }
    srand(time(NULL));
    for (int i = 0; i < len_data; i++) {
        time_t current_time = time(NULL);
        check = rand() / (double) RAND_MAX;
        for (int j = 0; j < n; j++) {
            if (check < laplas[j]) {
                data[i] += laplas_x[j];
                if (data[i] > quants[levels - 1]) {
                    data[i] = quants[levels - 1];
                }
                else if (data[i] < quants[0]) {
                    data[i] = quants[0];
                }
                break;
            }
        }
    }
}


void make_histogram1d(double* data, int len_data, double* hist_data, double* quants, int levels) {
    for (int i = 0; i < len_data; i++) {
        hist_data[(int)((data[i] - quants[0]) / (quants[1] - quants[0]))]++;
    }
}


void make_histogram2d(double** data, int len_data_x, int len_data_y, double* hist_data, double* quants, int levels) {
    for (int i = 0; i < len_data_y; i++) {
        for (int k = 0; k < len_data_x; k++) {
            for (int j = 0; j < levels; j++) {
                if (data[i][k] < quants[j]) {
                    hist_data[j]++;
                    break;
                }
            }
        }
    }
}


double mean(unsigned long int* hist_data, double* quants, int levels) {
    double sum = 0;
    double summa = 0;
    for (int i = 0; i < levels; i++) {
        sum += quants[i] * hist_data[i];
        summa += hist_data[i];
    }
    return sum / summa;
}


double variance(unsigned long int* hist_data, double* quants, int levels) {
    double sum = 0;
    double summa = 0;
    double me = mean(hist_data, quants, levels);
    for (int i = 0; i < levels; i++) {
        sum += pow(quants[i], 2) * hist_data[i];
        summa += hist_data[i];
    }
    return sum / summa - me * me;
}


double quartile_dist(unsigned long int* hist_data, double* quants, int levels) {
    double summa = 0;
    for (int i = 0; i < levels; i++) {
        summa += hist_data[i];
    }

    double tmp_summa = 0;
    int left_point, right_point;
    for (int i = 0; i < levels; i++) {
        tmp_summa += hist_data[i];
        if (tmp_summa > summa / 4.) {
            left_point = i;
            break;
        }
    }
    for (int i = left_point + 1; i < levels; i++) {
        tmp_summa += hist_data[i];
        if (tmp_summa > 3 * summa / 4.) {
            right_point = i;
        }
    }
    return quants[right_point] - quants[left_point];
}


void ImageProcessingGray(unsigned char* pOut
			 , unsigned char* pIn
			 , size_t nWidth
			 , size_t nHeight)
{
	unsigned long int histogram[256];
	for (int i = 0; i < 256; i++) {
		histogram[i] = 0;
	}

	for (size_t y = 0; y < nHeight; ++y)
	{
		for (size_t x = 0; x < nWidth; ++x)
		{

			histogram[pIn[nWidth * y + x]] += 1;
		}
	}
    double quants[256];
    for (int i = 0; i < 256; i++) {
        quants[i] = i;
    }

	printf("Bias = %f\n", mean(histogram, quants, 256));
	printf("Variance = %f\n", variance(histogram, quants, 256));
	printf("Quartile distance = %u\n", quartile_dist(histogram, quants, 256));

	return;
}


int main() {
    const int len_data = 1000;
    double data[len_data];
    double h = 0.01;
    for (int i = 0; i < len_data; i++) {
        data[i] = my_sin(i * h, 5, 0);
    }

    const int levels = 256;
    double right = INT64_MIN, left = INT64_MAX;
    for (int i = 0; i < len_data; i++) {
        if (data[i] > right) {
            right = data[i];
        }
        if (data[i] < left) {
            left = data[i];
        }
    }   

    double quants[levels];
    for (int i = 0; i < levels; i++) {
        quants[i] = left + i * (right - left) / (levels - 1); 
    }

    quantization(data, len_data, quants, levels);
    add_Gauss_noise(10, data, len_data, quants, levels);

    const char csv_file_name[64] = "work3.csv";
    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal\n";
    for (size_t i = 0; i < 1000; ++i){
        csv_file << (i * h) << "," << data[i] << "\n";
    }
    csv_file.close();
}