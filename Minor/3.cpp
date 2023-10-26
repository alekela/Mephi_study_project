#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <random>
#include "windows.h"


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


void add_impulse_noise(double P, double* data, int len_data) {
    double check;
    double a = rand();
    for (int i = 0; i < len_data; i++) {
        check = rand() / INT64_MAX;
        if (check < P) {
            data[i] = a;
        }
    }
}


double gauss(double x, double sigma) {
    return 1. / sqrt(2 * M_PI) / sigma * exp(- (x * x / 2. / sigma / sigma));
}


void add_Gauss_noise(double disp, double* data, int len_data) {
    double check;
    int n = 100;
    double laplas[n];
    double laplas_x[n];
    double h = 6. * disp / (n - 1);
    for (int i = 0; i< n; i++) {
        laplas_x[i] = -3. * disp + i * h;
        if (i == 0) {
            laplas[i] = 0;
        }
        else {
        laplas[i] = laplas[i - 1] + gauss(laplas_x[i], disp);
        }
    }
    for (int i = 0; i < len_data; i++) {
        check = rand() / INT64_MAX;
        for (int j = 0; j < n; j++) {
            if (check < laplas[j]) {
                data[i] += laplas_x[j];
            }
        }
    }
}


void make_histogram1d(double* data, int len_data, double* hist_data, double* quants, int levels) {
    for (int i = 0; i < len_data; i++) {
        for (int j = 0; j < levels; j++) {
            if (data[i] == quants[j]) {
                hist_data[j]++;
                break;
            }
        }
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


double mean(double* hist_data, double* quants, int levels) {
    double sum = 0;
    double summa = 0;
    for (int i = 0; i < levels; i++) {
        sum += quants[i] * hist_data[i];
        summa += hist_data[i];
    }
    return sum / summa;
}


double variance(double* hist_data, double* quants, int levels) {
    double sum = 0;
    double summa = 0;
    double me = mean(hist_data, quants, levels);
    for (int i = 0; i < levels; i++) {
        sum += pow(quants[i], 2) * hist_data[i];
        summa += hist_data[i];
    }
    return sum / summa - me * me;
}


double quartile_dist(double* hist_data, double* quants, int levels) {
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


int main() {
    int len_data = 1000;
    double data[len_data];
    double h = 0.01;
    for (int i = 0; i < len_data; i++) {
        data[i] = my_sin(i * h, 5, 0);
    }

    int levels = 256;
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

    add_impulse_noise(0.1, data, 1000);
    const char csv_file_name[64] = "work3.csv";
    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal\n";
    for (size_t i = 0; i < 1000; ++i){
        csv_file << (i * h) << "," << data[i] << "\n";
    }
    csv_file.close();
}
